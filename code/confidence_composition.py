"""Confidence-aware empirical adaptive composition.

This module implements the full propose-then-verify shape of Algorithm 1 plus
the adaptive composition loop from Algorithm 2.

Compared with ``empirical_composition.py``, this version uses two independent
sample supports:

* a training support to choose the Gaussian noise scale;
* a validation support to certify the centered local condition with a
  one-sided Hoeffding margin.

The code is still simulation-oriented: it uses a finite sampled support and a
Gaussian reference with an analytic alpha-divergence.  The validation check is
the operational version of the branch-local condition

    E[A_{t-1}(X) (D_t(X) - r_t)] <= 0.

If the empirical validation mean plus the Hoeffding margin is at most zero, the
candidate round mechanism is accepted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Sequence

import math
import numpy as np

from empirical_composition import (
    Array,
    LeakageFn,
    as_1d_output,
    as_2d_array,
    gaussian_logpdf,
    logsumexp,
    normalized_from_log_weights,
    safe_exp,
)


@dataclass
class ValidationStats:
    """One-sided Hoeffding validation statistics for a candidate."""

    mean: float
    lower: float
    upper: float
    beta: float
    upper_confidence_mean: float
    accepted: bool


@dataclass
class ConfidenceRoundResult:
    """Transcript and validation summary for one accepted round."""

    t: int
    query: Any
    sigma: float
    log_r_t: float
    r_t: float
    candidates_tried: int
    validation: ValidationStats
    reference_mean: Array
    raw_output: Array
    noisy_output: Array
    log_delta_bound: float
    delta_bound: float


def hoeffding_beta(
    lower: float,
    upper: float,
    sample_count: int,
    gamma: float,
    max_candidates: int,
) -> float:
    """One-sided validation margin with a union bound over candidates."""

    if sample_count <= 0:
        raise ValueError("sample_count must be positive")
    if not (0.0 < gamma < 1.0):
        raise ValueError("gamma must be in (0, 1)")
    if max_candidates <= 0:
        raise ValueError("max_candidates must be positive")
    return (upper - lower) * math.sqrt(math.log(max_candidates / gamma) / (2.0 * sample_count))


def required_validation_samples(
    range_width: float,
    margin: float,
    gamma: float,
    max_candidates: int,
) -> int:
    """Samples needed to make the Hoeffding beta no larger than ``margin``."""

    if range_width < 0:
        raise ValueError("range_width must be nonnegative")
    if margin <= 0:
        raise ValueError("margin must be positive")
    value = (range_width * range_width * math.log(max_candidates / gamma)) / (2.0 * margin * margin)
    return int(math.ceil(value))


def log10_required_validation_samples(
    log_range_width: float,
    log_margin: float,
    gamma: float,
    max_candidates: int,
) -> float:
    """Log10 of required validation samples, useful when counts are enormous."""

    return (
        2.0 * (log_range_width - log_margin)
        + math.log(math.log(max_candidates / gamma) / 2.0)
    ) / math.log(10.0)


class ConfidenceAdaptiveComposer:
    """Confidence-aware Gaussian adaptive composer."""

    def __init__(
        self,
        train_samples: Sequence[Any],
        validation_samples: Sequence[Any],
        leakage_fn: LeakageFn,
        alpha: float,
        log_r_factors: Sequence[float],
        train_prior_log_probs: Optional[Sequence[float]] = None,
        validation_prior_log_probs: Optional[Sequence[float]] = None,
        gamma: float = 0.05,
        max_candidates: int = 20,
        training_slack: float = 0.8,
        candidate_increment_fraction: float = 0.15,
        rng: Optional[np.random.Generator] = None,
        sigma_floor: float = 1e-9,
        sigma_initial_hi: float = 1.0,
        sigma_growth: float = 2.0,
        sigma_tol: float = 1e-6,
        max_sigma: float = 1e9,
        max_binary_steps: int = 80,
    ) -> None:
        if alpha <= 1:
            raise ValueError("alpha must be greater than 1")
        if len(train_samples) == 0 or len(validation_samples) == 0:
            raise ValueError("training and validation samples must be nonempty")
        if not (0.0 < training_slack < 1.0):
            raise ValueError("training_slack must be in (0, 1)")

        self.train_samples = list(train_samples)
        self.validation_samples = list(validation_samples)
        self.leakage_fn = leakage_fn
        self.alpha = float(alpha)
        self.log_r_factors = list(log_r_factors)
        self.gamma = float(gamma)
        self.max_candidates = int(max_candidates)
        self.training_slack = float(training_slack)
        self.candidate_increment_fraction = float(candidate_increment_fraction)
        self.rng = rng if rng is not None else np.random.default_rng()
        self.sigma_floor = float(sigma_floor)
        self.sigma_initial_hi = float(sigma_initial_hi)
        self.sigma_growth = float(sigma_growth)
        self.sigma_tol = float(sigma_tol)
        self.max_sigma = float(max_sigma)
        self.max_binary_steps = int(max_binary_steps)

        self.train_weights = np.full(len(self.train_samples), 1.0 / len(self.train_samples))
        self.validation_weights = np.full(
            len(self.validation_samples),
            1.0 / len(self.validation_samples),
        )
        self.train_log_weights = np.log(self.train_weights)
        self.validation_log_weights = np.log(self.validation_weights)

        if train_prior_log_probs is None:
            self.train_prior_log_probs = self.train_log_weights.copy()
        else:
            self.train_prior_log_probs = np.asarray(train_prior_log_probs, dtype=float)
            if self.train_prior_log_probs.shape != (len(self.train_samples),):
                raise ValueError("train_prior_log_probs must match train_samples")
        if validation_prior_log_probs is None:
            self.validation_prior_log_probs = self.validation_log_weights.copy()
        else:
            self.validation_prior_log_probs = np.asarray(validation_prior_log_probs, dtype=float)
            if self.validation_prior_log_probs.shape != (len(self.validation_samples),):
                raise ValueError("validation_prior_log_probs must match validation_samples")

        self.reset()

    def reset(self) -> None:
        self.transcript: list[ConfidenceRoundResult] = []
        self.train_log_prefix = (self.alpha - 1.0) * self.train_prior_log_probs.copy()
        self.validation_log_prefix = (self.alpha - 1.0) * self.validation_prior_log_probs.copy()
        self._normalize_prefixes()
        self.log_product_r = 0.0

    def run(self, queries: Sequence[Any], actual_secret: Any) -> list[ConfidenceRoundResult]:
        if len(queries) != len(self.log_r_factors):
            raise ValueError("queries must have one entry per round")

        for t, (query, log_r_t) in enumerate(zip(queries, self.log_r_factors), start=1):
            train_means = self._evaluate(self.train_samples, query, t)
            validation_means = self._evaluate(self.validation_samples, query, t)
            sigma, reference_mean, stats, candidates = self.calibrate_round(
                train_means,
                validation_means,
                log_r_t,
            )

            raw_output = as_1d_output(self.leakage_fn(actual_secret, query, t, self.transcript))
            noise = self.rng.normal(0.0, sigma, size=raw_output.shape)
            noisy_output = raw_output + noise

            self._update_prefix(
                train_means,
                validation_means,
                reference_mean,
                sigma,
                noisy_output,
            )
            self.log_product_r += log_r_t

            result = ConfidenceRoundResult(
                t=t,
                query=query,
                sigma=sigma,
                log_r_t=log_r_t,
                r_t=safe_exp(log_r_t),
                candidates_tried=candidates,
                validation=stats,
                reference_mean=reference_mean.copy(),
                raw_output=raw_output,
                noisy_output=noisy_output,
                log_delta_bound=self._log_delta_bound(),
                delta_bound=safe_exp(self._log_delta_bound()),
            )
            self.transcript.append(result)

        return self.transcript

    def calibrate_round(
        self,
        train_means: Array,
        validation_means: Array,
        log_r_t: float,
    ) -> tuple[float, Array, ValidationStats, int]:
        """Propose on training samples and verify on validation samples."""

        if log_r_t <= 0:
            raise ValueError("finite-noise validation requires log_r_t > 0")

        reference_mean = self._reference_mean(train_means, self.train_log_prefix, self.train_log_weights)
        train_log_budget = log_r_t + math.log(self.training_slack)
        sigma = self._binary_search_sigma(train_means, reference_mean, train_log_budget)

        last_stats: Optional[ValidationStats] = None
        for candidate_index in range(1, self.max_candidates + 1):
            stats = self._validation_stats(validation_means, reference_mean, sigma, log_r_t)
            last_stats = stats
            if stats.accepted:
                return sigma, reference_mean, stats, candidate_index

            if candidate_index < self.max_candidates:
                increment = max(self.sigma_floor, self.candidate_increment_fraction * sigma)
                sigma = math.sqrt(sigma * sigma + increment * increment)

        if last_stats is None:
            raise RuntimeError("validation loop did not run")
        raise RuntimeError(
            "no candidate passed validation; increase training_slack, max_candidates, "
            "or validation sample size"
        )

    def _evaluate(self, samples: Sequence[Any], query: Any, t: int) -> Array:
        return as_2d_array(self.leakage_fn(sample, query, t, self.transcript) for sample in samples)

    def _reference_mean(self, means: Array, log_prefix: Array, log_weights: Array) -> Array:
        normalized = normalized_from_log_weights(log_weights + log_prefix)
        return normalized @ means

    def _binary_search_sigma(self, means: Array, reference_mean: Array, log_budget: float) -> float:
        sq_dist = np.sum((means - reference_mean.reshape(1, -1)) ** 2, axis=1)
        if float(np.max(sq_dist)) == 0.0:
            return self.sigma_floor

        lo = self.sigma_floor
        hi = max(self.sigma_initial_hi, lo)
        while self._log_local_ratio(sq_dist, hi, self.train_log_prefix, self.train_log_weights) > log_budget:
            hi *= self.sigma_growth
            if hi > self.max_sigma:
                raise RuntimeError("failed to find feasible sigma")

        for _ in range(self.max_binary_steps):
            mid = 0.5 * (lo + hi)
            if self._log_local_ratio(sq_dist, mid, self.train_log_prefix, self.train_log_weights) <= log_budget:
                hi = mid
            else:
                lo = mid
            if hi - lo <= self.sigma_tol * max(1.0, hi):
                break
        return hi

    def _log_local_ratio(
        self,
        sq_dist: Array,
        sigma: float,
        log_prefix: Array,
        log_weights: Array,
    ) -> float:
        log_div = self.alpha * (self.alpha - 1.0) * sq_dist / (2.0 * sigma * sigma)
        log_num = logsumexp(log_weights + log_prefix + log_div)
        log_den = logsumexp(log_weights + log_prefix)
        return log_num - log_den

    def _validation_stats(
        self,
        means: Array,
        reference_mean: Array,
        sigma: float,
        log_r_t: float,
    ) -> ValidationStats:
        sq_dist = np.sum((means - reference_mean.reshape(1, -1)) ** 2, axis=1)
        log_div = self.alpha * (self.alpha - 1.0) * sq_dist / (2.0 * sigma * sigma)
        if float(np.max(log_div)) > 700 or log_r_t > 700:
            raise OverflowError(
                "validation values exceed float range; use log-scale sample-size "
                "diagnostics or a smaller alpha"
            )

        prefix = np.exp(self.validation_log_prefix)
        div = np.exp(log_div)
        r_t = math.exp(log_r_t)
        values = prefix * (div - r_t)
        mean = float(np.mean(values))
        lower = float(np.min(values))
        upper = float(np.max(values))
        beta = hoeffding_beta(
            lower,
            upper,
            len(values),
            self.gamma,
            self.max_candidates,
        )
        upper_confidence_mean = mean + beta
        return ValidationStats(
            mean=mean,
            lower=lower,
            upper=upper,
            beta=beta,
            upper_confidence_mean=upper_confidence_mean,
            accepted=upper_confidence_mean <= 0.0,
        )

    def _update_prefix(
        self,
        train_means: Array,
        validation_means: Array,
        reference_mean: Array,
        sigma: float,
        noisy_output: Array,
    ) -> None:
        logp_reference = gaussian_logpdf(noisy_output, reference_mean, sigma)[0]
        self.train_log_prefix = self.train_log_prefix + self.alpha * (
            gaussian_logpdf(noisy_output, train_means, sigma) - logp_reference
        )
        self.validation_log_prefix = self.validation_log_prefix + self.alpha * (
            gaussian_logpdf(noisy_output, validation_means, sigma) - logp_reference
        )
        self._normalize_prefixes()

    def _normalize_prefixes(self) -> None:
        self.train_log_prefix -= logsumexp(self.train_log_weights + self.train_log_prefix)
        self.validation_log_prefix -= logsumexp(
            self.validation_log_weights + self.validation_log_prefix
        )

    def _log_delta_bound(self) -> float:
        log_base = logsumexp(
            self.train_log_weights + (self.alpha - 1.0) * self.train_prior_log_probs
        )
        return (log_base + self.log_product_r) / self.alpha
