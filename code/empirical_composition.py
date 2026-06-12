"""Empirical black-box privatization with adaptive composition accounting.

This module implements a first simulation version of Algorithm 2 in the paper:
"Randomization Optimizer under Adversarial Composition".

Simplification used here:
    A sampled support is treated as the whole secret distribution.  The prior is
    the empirical distribution on those samples, so expectations over D are exact
    finite sums.  We deliberately omit Hoeffding/confidence parameters and the
    independent verification step from Algorithm 1.

The implementation uses additive isotropic Gaussian output noise.  At each
round, it picks the smallest sigma found by binary search such that the empirical
branch-local alpha-divergence condition holds:

    E_emp[A_{t-1}(X) D_alpha(M_t(X) || W_t)] <= r_t E_emp[A_{t-1}(X)].

For this first simulation, W_t is a Gaussian with the same covariance as the
mechanism and center equal to the prefix-weighted empirical leakage mean.  This
choice keeps D_alpha between the Gaussian conditionals analytic.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, List, Optional, Protocol, Sequence

import math
import numpy as np


Array = np.ndarray
LeakageFn = Callable[[Any, Any, int, Sequence["RoundResult"]], Any]


class Adversary(Protocol):
    """Minimal adversary interface used by the adaptive loop."""

    def next_query(self, transcript: Sequence["RoundResult"]) -> Any:
        """Return the next query after observing the current transcript."""


@dataclass
class RoundResult:
    """Public transcript and accounting summary for one adaptive round."""

    t: int
    query: Any
    sigma: float
    r_t: float
    log_r_t: float
    reference_mean: Array
    raw_output: Array
    noisy_output: Array
    noise: Array
    log_local_ratio: float
    local_ratio: float
    log_delta_bound: float
    delta_bound: float


@dataclass
class CalibrationResult:
    """The result of one empirical local calibration problem."""

    sigma: float
    reference_mean: Array
    log_local_ratio: float
    local_ratio: float
    max_sq_distance: float


def as_2d_array(values: Iterable[Any]) -> Array:
    """Convert scalar/vector leakage evaluations to an (m, d) float array."""

    rows = [np.asarray(v, dtype=float).reshape(-1) for v in values]
    if not rows:
        raise ValueError("at least one sample is required")
    width = rows[0].shape[0]
    if any(row.shape[0] != width for row in rows):
        raise ValueError("all leakage outputs must have the same dimension")
    return np.vstack(rows)


def as_1d_output(value: Any) -> Array:
    """Convert one leakage output to a flat float vector."""

    return np.asarray(value, dtype=float).reshape(-1)


def logsumexp(log_values: Array) -> float:
    """Small local logsumexp to avoid depending on scipy."""

    values = np.asarray(log_values, dtype=float)
    max_value = float(np.max(values))
    if not math.isfinite(max_value):
        return max_value
    return max_value + math.log(float(np.sum(np.exp(values - max_value))))


def normalized_from_log_weights(log_weights: Array) -> Array:
    """Return normalized positive weights represented in log-space."""

    log_total = logsumexp(log_weights)
    return np.exp(log_weights - log_total)


def safe_exp(value: float) -> float:
    """Exponentiate for reporting without raising on very large values."""

    if value > 700:
        return float("inf")
    if value < -745:
        return 0.0
    return math.exp(value)


def gaussian_logpdf(y: Array, means: Array, sigma: float) -> Array:
    """Log-density of N(mean, sigma^2 I) at y for each row in means."""

    if sigma <= 0:
        raise ValueError("sigma must be positive")
    y_vec = np.asarray(y, dtype=float).reshape(1, -1)
    mu = np.asarray(means, dtype=float)
    if mu.ndim == 1:
        mu = mu.reshape(1, -1)
    if mu.shape[1] != y_vec.shape[1]:
        raise ValueError("mean and observation dimensions do not match")
    dim = mu.shape[1]
    sq_norm = np.sum((mu - y_vec) ** 2, axis=1)
    log_norm = -0.5 * dim * math.log(2.0 * math.pi * sigma * sigma)
    return log_norm - sq_norm / (2.0 * sigma * sigma)


def empirical_base_log(
    weights: Array,
    alpha: float,
    prior_log_probs: Optional[Array] = None,
) -> float:
    """Compute log E_D[pi(X)^(alpha-1)] on the sampled support.

    ``weights`` are the sampling/expectation weights over the finite support.
    ``prior_log_probs`` are the log prior probabilities pi(x) that appear in the
    PRW bound.  By default they are the same distribution, which recovers the
    pure empirical-prior simplification.
    """

    log_weights = np.log(np.asarray(weights, dtype=float))
    if prior_log_probs is None:
        prior_log_probs = log_weights
    return logsumexp(log_weights + (alpha - 1.0) * np.asarray(prior_log_probs, dtype=float))


def equal_r_factors_for_target(
    weights: Array,
    alpha: float,
    target_delta: float,
    rounds: int,
    prior_log_probs: Optional[Array] = None,
) -> List[float]:
    """Use equal per-round r_t factors for a desired PRW alpha bound.

    The composition theorem gives

        delta^alpha <= E_D[pi(X)^(alpha-1)] prod_t r_t.

    With equal factors, each r_t is

        (target_delta^alpha / base)^(1 / rounds).
    """

    if rounds <= 0:
        raise ValueError("rounds must be positive")
    if not (0.0 < target_delta <= 1.0):
        raise ValueError("target_delta must be in (0, 1]")
    log_base = empirical_base_log(weights, alpha, prior_log_probs=prior_log_probs)
    log_r = (alpha * math.log(target_delta) - log_base) / rounds
    r = math.exp(log_r)
    if r < 1.0:
        raise ValueError(
            "target_delta is below the no-leakage empirical alpha bound; "
            f"equal factor would be {r:.6g} < 1"
        )
    return [r for _ in range(rounds)]


class EmpiricalCompositionAccountant:
    """Tracks the theorem-level composition bound."""

    def __init__(
        self,
        weights: Array,
        alpha: float,
        prior_log_probs: Optional[Array] = None,
    ) -> None:
        self.weights = np.asarray(weights, dtype=float)
        self.alpha = float(alpha)
        if prior_log_probs is None:
            self.prior_log_probs = np.log(self.weights)
        else:
            self.prior_log_probs = np.asarray(prior_log_probs, dtype=float)
        self.log_base = empirical_base_log(
            self.weights,
            self.alpha,
            prior_log_probs=self.prior_log_probs,
        )
        self.log_product_r = 0.0

    def add_round(self, r_t: float) -> None:
        if r_t <= 0:
            raise ValueError("r_t must be positive")
        self.add_log_round(math.log(r_t))

    def add_log_round(self, log_r_t: float) -> None:
        if not math.isfinite(log_r_t):
            raise ValueError("log_r_t must be finite")
        self.log_product_r += log_r_t

    @property
    def log_delta_bound(self) -> float:
        return (self.log_base + self.log_product_r) / self.alpha

    @property
    def delta_bound(self) -> float:
        return safe_exp(self.log_delta_bound)


class EmpiricalAdaptiveComposer:
    """Run empirical adaptive Gaussian calibration and composition accounting."""

    def __init__(
        self,
        samples: Sequence[Any],
        leakage_fn: LeakageFn,
        alpha: float,
        r_factors: Optional[Sequence[float]] = None,
        log_r_factors: Optional[Sequence[float]] = None,
        weights: Optional[Sequence[float]] = None,
        prior_log_probs: Optional[Sequence[float]] = None,
        rng: Optional[np.random.Generator] = None,
        sigma_floor: float = 1e-9,
        sigma_initial_hi: float = 1.0,
        sigma_growth: float = 2.0,
        sigma_tol: float = 1e-6,
        max_sigma: float = 1e9,
        max_binary_steps: int = 80,
        normalize_prefix: bool = True,
    ) -> None:
        if alpha <= 1:
            raise ValueError("alpha must be greater than 1")
        if len(samples) == 0:
            raise ValueError("samples must be nonempty")
        if r_factors is not None and log_r_factors is not None:
            raise ValueError("provide either r_factors or log_r_factors, not both")
        self.samples = list(samples)
        self.leakage_fn = leakage_fn
        self.alpha = float(alpha)
        self.r_factors = list(r_factors) if r_factors is not None else None
        self.log_r_factors = list(log_r_factors) if log_r_factors is not None else None
        self.rng = rng if rng is not None else np.random.default_rng()
        self.sigma_floor = float(sigma_floor)
        self.sigma_initial_hi = float(sigma_initial_hi)
        self.sigma_growth = float(sigma_growth)
        self.sigma_tol = float(sigma_tol)
        self.max_sigma = float(max_sigma)
        self.max_binary_steps = int(max_binary_steps)
        self.normalize_prefix = bool(normalize_prefix)

        if weights is None:
            self.weights = np.full(len(self.samples), 1.0 / len(self.samples))
        else:
            self.weights = np.asarray(weights, dtype=float)
            if self.weights.shape != (len(self.samples),):
                raise ValueError("weights must have one entry per sample")
            if np.any(self.weights <= 0):
                raise ValueError("weights must be strictly positive")
            self.weights = self.weights / float(np.sum(self.weights))

        self.log_weights = np.log(self.weights)
        if prior_log_probs is None:
            self.prior_log_probs = self.log_weights.copy()
        else:
            self.prior_log_probs = np.asarray(prior_log_probs, dtype=float)
            if self.prior_log_probs.shape != (len(self.samples),):
                raise ValueError("prior_log_probs must have one entry per sample")
        self.reset()

    def reset(self) -> None:
        """Reset prefix coefficients and theorem-level accountant."""

        self.transcript: List[RoundResult] = []
        self.log_prefix = (self.alpha - 1.0) * self.prior_log_probs.copy()
        self.accountant = EmpiricalCompositionAccountant(
            self.weights,
            self.alpha,
            prior_log_probs=self.prior_log_probs,
        )

    def run(
        self,
        adversary: Adversary | Callable[[Sequence[RoundResult]], Any],
        rounds: Optional[int] = None,
        actual_secret: Optional[Any] = None,
    ) -> List[RoundResult]:
        """Run the empirical adaptive composition mechanism."""

        if self.r_factors is None and rounds is None:
            if self.log_r_factors is None:
                raise ValueError("rounds is required when r_factors/log_r_factors is not set")
        if self.r_factors is not None:
            total_rounds = len(self.r_factors)
        elif self.log_r_factors is not None:
            total_rounds = len(self.log_r_factors)
        else:
            total_rounds = int(rounds)
        if self.r_factors is not None and rounds is not None and rounds != len(self.r_factors):
            raise ValueError("rounds must match len(r_factors)")
        if self.log_r_factors is not None and rounds is not None and rounds != len(self.log_r_factors):
            raise ValueError("rounds must match len(log_r_factors)")

        if actual_secret is None:
            actual_secret = self.samples[int(self.rng.integers(len(self.samples)))]

        for t in range(1, total_rounds + 1):
            if self.r_factors is not None:
                r_t = self.r_factors[t - 1]
                log_r_t = math.log(r_t)
            elif self.log_r_factors is not None:
                log_r_t = self.log_r_factors[t - 1]
                r_t = safe_exp(log_r_t)
            else:
                log_r_t = 0.0
                r_t = 1.0
            query = self._next_query(adversary)
            means = self._evaluate_sample_means(query, t)
            calibration = self.calibrate_round(means, log_r_t=log_r_t)

            raw_output = as_1d_output(self.leakage_fn(actual_secret, query, t, self.transcript))
            if raw_output.shape[0] != means.shape[1]:
                raise ValueError("actual leakage output dimension differs from sample outputs")
            noise = self.rng.normal(0.0, calibration.sigma, size=raw_output.shape)
            noisy_output = raw_output + noise

            self._update_prefix(means, calibration.reference_mean, calibration.sigma, noisy_output)
            self.accountant.add_log_round(log_r_t)

            result = RoundResult(
                t=t,
                query=query,
                sigma=calibration.sigma,
                r_t=r_t,
                log_r_t=log_r_t,
                reference_mean=calibration.reference_mean.copy(),
                raw_output=raw_output,
                noisy_output=noisy_output,
                noise=noise,
                log_local_ratio=calibration.log_local_ratio,
                local_ratio=calibration.local_ratio,
                log_delta_bound=self.accountant.log_delta_bound,
                delta_bound=self.accountant.delta_bound,
            )
            self.transcript.append(result)

        return self.transcript

    def calibrate_round(
        self,
        means: Array,
        r_t: Optional[float] = None,
        log_r_t: Optional[float] = None,
    ) -> CalibrationResult:
        """Find the smallest Gaussian sigma satisfying the local condition."""

        if log_r_t is None:
            if r_t is None:
                raise ValueError("either r_t or log_r_t is required")
            if r_t <= 0:
                raise ValueError("r_t must be positive")
            log_r_t = math.log(r_t)
        if log_r_t < 0.0:
            raise ValueError("r_t must be at least 1 for Gaussian smoothing to be feasible")
        reference_mean = self._reference_mean(means)
        sq_dist = np.sum((means - reference_mean.reshape(1, -1)) ** 2, axis=1)
        max_sq_distance = float(np.max(sq_dist))

        if max_sq_distance == 0.0:
            sigma = self.sigma_floor
            log_ratio = self._log_local_ratio_from_sqdist(sq_dist, sigma)
            return CalibrationResult(
                sigma=sigma,
                reference_mean=reference_mean,
                log_local_ratio=log_ratio,
                local_ratio=safe_exp(log_ratio),
                max_sq_distance=max_sq_distance,
            )

        if log_r_t == 0.0:
            raise ValueError("r_t=1 is feasible only when all conditional means match")

        lo = self.sigma_floor
        hi = max(self.sigma_initial_hi, lo)
        while not self._passes_from_sqdist(sq_dist, hi, log_r_t):
            hi *= self.sigma_growth
            if hi > self.max_sigma:
                raise RuntimeError(
                    "failed to find a feasible sigma before max_sigma; "
                    "try larger r_t or max_sigma"
                )

        for _ in range(self.max_binary_steps):
            mid = 0.5 * (lo + hi)
            if self._passes_from_sqdist(sq_dist, mid, log_r_t):
                hi = mid
            else:
                lo = mid
            if hi - lo <= self.sigma_tol * max(1.0, hi):
                break

        sigma = hi
        log_ratio = self._log_local_ratio_from_sqdist(sq_dist, sigma)
        return CalibrationResult(
            sigma=sigma,
            reference_mean=reference_mean,
            log_local_ratio=log_ratio,
            local_ratio=safe_exp(log_ratio),
            max_sq_distance=max_sq_distance,
        )

    def _next_query(
        self,
        adversary: Adversary | Callable[[Sequence[RoundResult]], Any],
    ) -> Any:
        if hasattr(adversary, "next_query"):
            return adversary.next_query(self.transcript)  # type: ignore[union-attr]
        return adversary(self.transcript)  # type: ignore[misc]

    def _evaluate_sample_means(self, query: Any, t: int) -> Array:
        values = [self.leakage_fn(sample, query, t, self.transcript) for sample in self.samples]
        return as_2d_array(values)

    def _reference_mean(self, means: Array) -> Array:
        """Prefix-weighted empirical center for the Gaussian reference."""

        log_importance = self.log_weights + self.log_prefix
        normalized = normalized_from_log_weights(log_importance)
        return normalized @ means

    def _log_local_ratio_from_sqdist(self, sq_dist: Array, sigma: float) -> float:
        """Log of E[A D_alpha] / E[A] for equal-covariance Gaussians."""

        if sigma <= 0:
            return float("inf")
        log_div = self.alpha * (self.alpha - 1.0) * sq_dist / (2.0 * sigma * sigma)
        log_importance = self.log_weights + self.log_prefix
        log_num = logsumexp(log_importance + log_div)
        log_den = logsumexp(log_importance)
        return log_num - log_den

    def _passes_from_sqdist(self, sq_dist: Array, sigma: float, log_r_t: float) -> bool:
        return self._log_local_ratio_from_sqdist(sq_dist, sigma) <= log_r_t

    def _update_prefix(
        self,
        means: Array,
        reference_mean: Array,
        sigma: float,
        noisy_output: Array,
    ) -> None:
        """Update A_t up to a common positive scale."""

        logp_candidates = gaussian_logpdf(noisy_output, means, sigma)
        logp_reference = gaussian_logpdf(noisy_output, reference_mean, sigma)[0]
        self.log_prefix = self.log_prefix + self.alpha * (logp_candidates - logp_reference)

        if self.normalize_prefix:
            # Scaling all A_t values by a common constant does not change any
            # future local condition.  This keeps the log-prefix numerically sane.
            self.log_prefix = self.log_prefix - logsumexp(self.log_weights + self.log_prefix)


class FixedQueryAdversary:
    """A simple adversary useful for smoke tests and scripted experiments."""

    def __init__(self, queries: Sequence[Any]) -> None:
        if len(queries) == 0:
            raise ValueError("queries must be nonempty")
        self.queries = list(queries)

    def next_query(self, transcript: Sequence[RoundResult]) -> Any:
        index = len(transcript)
        if index >= len(self.queries):
            raise IndexError("no query configured for this round")
        return self.queries[index]


def toy_linear_leakage(secret: Any, query: Any, t: int, transcript: Sequence[RoundResult]) -> Array:
    """Toy scalar leakage: normalized dot product between a bit vector and query."""

    secret_vec = np.asarray(secret, dtype=float).reshape(-1)
    query_vec = np.asarray(query, dtype=float).reshape(-1)
    if secret_vec.shape != query_vec.shape:
        raise ValueError("secret and query vectors must have the same shape")
    return np.array([float(np.dot(secret_vec, query_vec) / math.sqrt(secret_vec.size))])


def demo() -> None:
    """Run a tiny empirical composition simulation."""

    rng = np.random.default_rng(7)
    sample_count = 256
    secret_bits = 16
    rounds = 5
    alpha = 8.0
    samples = rng.integers(0, 2, size=(sample_count, secret_bits))
    queries = rng.choice([-1.0, 1.0], size=(rounds, secret_bits))
    weights = np.full(sample_count, 1.0 / sample_count)
    r_factors = equal_r_factors_for_target(
        weights=weights,
        alpha=alpha,
        target_delta=0.05,
        rounds=rounds,
    )

    composer = EmpiricalAdaptiveComposer(
        samples=samples,
        leakage_fn=toy_linear_leakage,
        alpha=alpha,
        r_factors=r_factors,
        rng=rng,
    )
    transcript = composer.run(FixedQueryAdversary(queries), actual_secret=samples[0])

    print("round  sigma       local_ratio   r_t          delta_bound")
    for row in transcript:
        print(
            f"{row.t:5d}  "
            f"{row.sigma:10.6f}  "
            f"{row.local_ratio:11.6f}  "
            f"{row.r_t:11.6f}  "
            f"{row.delta_bound:11.6f}"
        )


if __name__ == "__main__":
    demo()
