"""Empirical RSA timing figures with no closed-form ground truth.

The RSA timing channel is treated as a black-box simulator.  We sample exponent
strings, evaluate the timing proxy, discretize the scalar timing output into
integer ns bins, and compute empirical PRW/PAC quantities from those samples.
No Hoeffding/confidence correction and no closed-form ground-truth curve are
used in this first simulation pass.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import cvxpy as cp
import numpy as np
from PIL import Image, ImageDraw
from scipy.optimize import linprog
from scipy.special import gammaln, logsumexp
from scipy.sparse import lil_matrix

from empirical_composition import as_2d_array
from gaussian_comparison_figures import (
    font,
    invert_alpha_bound_logrho_vector,
    invert_kl_bound,
)
from rsa_timing_one_round import make_query, rsa_timing_proxy


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "rsa_timing_empirical"
SAMPLE_COUNT = 100_000
SEED = 20260606
PRW_ALPHAS = (16.0, 64.0, 256.0)
OPT_ALPHA = 64.0
NOISE_RADIUS = 2048
COMPOSITION_TRIALS = 5
COMPOSITION_MAX_ROUNDS = 50
COMPOSITION_ROUNDS = list(range(1, COMPOSITION_MAX_ROUNDS + 1))
COMPOSITION_X_TICKS = [1, 10, 20, 30, 40, 50]
COMPOSITION_PRW_ALPHA = 256.0
BOUNDED_LABEL = "Optimal zero-mean"
ONESIDED_LABEL = "Optimal non-negative noise"


@dataclass
class TimingCase:
    label: str
    bits: int
    samples: np.ndarray
    leakage_raw: np.ndarray
    leakage_bins: np.ndarray
    leakage_values: np.ndarray
    leakage_counts: np.ndarray
    log_prior: np.ndarray
    query: dict[str, int]
    prior_success_log2: float
    task_log2_bonus: float = 0.0


def log_binom_counts(bits: int) -> np.ndarray:
    k = np.arange(bits + 1, dtype=float)
    return gammaln(bits + 1) - gammaln(k + 1) - gammaln(bits - k + 1)


def log2_hamming_ball_size(bits: int, radius: int) -> float:
    return log2(float(logsumexp(log_binom_counts(bits)[: radius + 1])))


def log2(x: float) -> float:
    return x / math.log(2.0)


def make_timing_case(
    bits: int,
    *,
    label: str,
    seed: int,
    bit_prob_one: float = 0.5,
    bin_width: float = 1.0,
) -> TimingCase:
    rng = np.random.default_rng(seed)
    if bit_prob_one == 0.5:
        samples = rng.integers(0, 2, size=(SAMPLE_COUNT, bits), dtype=np.uint8)
    else:
        samples = (rng.random((SAMPLE_COUNT, bits)) < bit_prob_one).astype(np.uint8)
    query = make_query(rng)
    raw = as_2d_array(rsa_timing_proxy(sample, query, 1, []) for sample in samples)[:, 0]
    bins = np.rint(raw / bin_width).astype(int)
    weights = samples.sum(axis=1)
    if bit_prob_one == 0.5:
        log_prior = np.full(SAMPLE_COUNT, -bits * math.log(2.0))
        prior_success_log2 = -float(bits)
    else:
        p = bit_prob_one
        log_prior = weights * math.log(p) + (bits - weights) * math.log1p(-p)
        prior_success_log2 = bits * math.log2(max(p, 1.0 - p))
    values, counts = np.unique(bins, return_counts=True)
    return TimingCase(label, bits, samples, raw, bins, values.astype(int), counts.astype(float), log_prior, query, prior_success_log2)


def make_timing_case_from_query(
    samples: np.ndarray,
    *,
    bits: int,
    label: str,
    query: dict[str, int],
    log_prior: np.ndarray,
    prior_success_log2: float,
    bin_width: float = 1.0,
) -> TimingCase:
    raw = as_2d_array(rsa_timing_proxy(sample, query, 1, []) for sample in samples)[:, 0]
    bins = np.rint(raw / bin_width).astype(int)
    values, counts = np.unique(bins, return_counts=True)
    return TimingCase(
        label=label,
        bits=bits,
        samples=samples,
        leakage_raw=raw,
        leakage_bins=bins,
        leakage_values=values.astype(int),
        leakage_counts=counts.astype(float),
        log_prior=log_prior,
        query=query,
        prior_success_log2=prior_success_log2,
    )


def make_fast_composition_query(rng: np.random.Generator, bits: int) -> dict[str, object]:
    """A vectorized RSA-style query for composition panels.

    The mask/weights stand in for public-input-dependent operand effects during
    private-key modular exponentiation.  The dominant term is still the
    square-and-multiply dependence on the secret exponent bits.
    """

    return {
        "mask": rng.integers(0, 2, size=bits, dtype=np.uint8),
        "weights": rng.normal(0.0, 1.0, size=bits).astype(float),
        "base_id": int(rng.integers(2, 1 << 63)),
    }


def rsa_timing_composition_raw(samples: np.ndarray, query: dict[str, object]) -> np.ndarray:
    bits = samples.shape[1]
    hw = samples.sum(axis=1).astype(float)
    mask = np.asarray(query["mask"], dtype=np.uint8)
    xor_hw = np.bitwise_xor(samples, mask).sum(axis=1).astype(float)
    weights = np.asarray(query["weights"], dtype=float)
    operand_term = samples.astype(float) @ weights / math.sqrt(float(bits))
    return (
        1.05 * float(bits)
        + 1.20 * hw
        + 0.04 * xor_hw
        + 0.15 * operand_term
    )


def make_fast_composition_case(
    samples: np.ndarray,
    *,
    bits: int,
    label: str,
    query: dict[str, object],
    log_prior: np.ndarray,
    prior_success_log2: float,
    bin_width: float = 1.0,
) -> TimingCase:
    raw = rsa_timing_composition_raw(samples, query)
    bins = np.rint(raw / bin_width).astype(int)
    values, counts = np.unique(bins, return_counts=True)
    return TimingCase(
        label=label,
        bits=bits,
        samples=samples,
        leakage_raw=raw,
        leakage_bins=bins,
        leakage_values=values.astype(int),
        leakage_counts=counts.astype(float),
        log_prior=log_prior,
        query={"base": int(query["base_id"])},
        prior_success_log2=prior_success_log2,
    )


def make_poisson_hw_timing_case(
    bits: int,
    *,
    label: str,
    seed: int,
    lam: float,
    bin_width: float = 1.0,
) -> TimingCase:
    rng = np.random.default_rng(seed)
    weights = np.arange(bits + 1, dtype=int)
    log_counts = log_binom_counts(bits)
    log_h_mass = weights.astype(float) * math.log(lam) - gammaln(weights.astype(float) + 1.0)
    log_h_mass = log_h_mass - float(logsumexp(log_h_mass))
    h_probs = np.exp(log_h_mass)
    sampled_weights = rng.choice(weights, size=SAMPLE_COUNT, p=h_probs)

    samples = np.zeros((SAMPLE_COUNT, bits), dtype=np.uint8)
    for weight in np.unique(sampled_weights):
        row_idx = np.flatnonzero(sampled_weights == weight)
        if weight == 0:
            continue
        if weight == bits:
            samples[row_idx, :] = 1
            continue
        random_scores = rng.random((row_idx.size, bits))
        chosen = np.argpartition(random_scores, int(weight) - 1, axis=1)[:, : int(weight)]
        samples[row_idx[:, None], chosen] = 1

    query = make_query(rng)
    raw = as_2d_array(rsa_timing_proxy(sample, query, 1, []) for sample in samples)[:, 0]
    bins = np.rint(raw / bin_width).astype(int)
    log_prior = log_h_mass[sampled_weights] - log_counts[sampled_weights]
    prior_success_log2 = log2(float(np.max(log_h_mass - log_counts)))
    values, counts = np.unique(bins, return_counts=True)
    return TimingCase(label, bits, samples, raw, bins, values.astype(int), counts.astype(float), log_prior, query, prior_success_log2)


def empirical_log_weights(case: TimingCase) -> np.ndarray:
    return np.full(case.leakage_bins.size, -math.log(case.leakage_bins.size))


def grouped_probabilities(case: TimingCase) -> Tuple[np.ndarray, np.ndarray]:
    probs = case.leakage_counts.astype(float) / float(np.sum(case.leakage_counts))
    return case.leakage_values.astype(float), probs


def grouped_prw_weights(case: TimingCase, alpha: float) -> Tuple[np.ndarray, np.ndarray, float]:
    values = case.leakage_values
    logw = empirical_log_weights(case)
    if float(np.max(case.log_prior) - np.min(case.log_prior)) <= 1e-12:
        probs = case.leakage_counts.astype(float) / float(np.sum(case.leakage_counts))
        log_base = float((alpha - 1.0) * case.log_prior[0])
        return values.astype(float), probs, log_base

    inverse = np.searchsorted(values, case.leakage_bins)
    logw = empirical_log_weights(case)
    logs = []
    for idx in range(values.size):
        mask = inverse == idx
        logs.append(float(logsumexp(logw[mask] + (alpha - 1.0) * case.log_prior[mask])))
    log_group = np.asarray(logs, dtype=float)
    log_base = float(logsumexp(logw + (alpha - 1.0) * case.log_prior))
    normalized = np.exp(log_group - log_base)
    normalized /= float(np.sum(normalized))
    return values.astype(float), normalized, log_base


def integration_grid(values: np.ndarray, sigma: float) -> Tuple[np.ndarray, float]:
    lo = float(np.min(values)) - 8.0 * sigma
    hi = float(np.max(values)) + 8.0 * sigma
    dx = max(0.2, min(1.0, 0.04 * sigma))
    count = int(math.ceil((hi - lo) / dx)) + 1
    y = np.linspace(lo, hi, count)
    return y, float((hi - lo) / (count - 1))


def gaussian_logpdf_grid(y: np.ndarray, means: np.ndarray, sigma: float) -> np.ndarray:
    return -0.5 * math.log(2.0 * math.pi * sigma * sigma) - (
        (y[:, None] - means[None, :]) ** 2
    ) / (2.0 * sigma * sigma)


def mutual_information(case: TimingCase, sigma: float) -> float:
    values, probs = grouped_probabilities(case)
    y, dx = integration_grid(values, sigma)
    logpdf = gaussian_logpdf_grid(y, values, sigma)
    log_joint = np.log(probs)[None, :] + logpdf
    log_py = logsumexp(log_joint, axis=1)
    integrand = np.sum(np.exp(log_joint) * (logpdf - log_py[:, None]), axis=1)
    return float(np.sum(integrand) * dx)


def pac_alpha_best_log2(case: TimingCase, sigma: float, alpha_grid: np.ndarray) -> Tuple[float, float]:
    values, probs = grouped_probabilities(case)
    mu = float(np.sum(values * probs))
    q_task = 2.0 ** case.prior_success_log2
    best_log2 = 0.0
    best_alpha = float(alpha_grid[0])

    def evaluate(alpha_values: np.ndarray) -> None:
        nonlocal best_log2, best_alpha
        alpha_values = np.unique(alpha_values[alpha_values > 1.0])
        if alpha_values.size == 0:
            return
        c = alpha_values * (alpha_values - 1.0) / 2.0
        scaled_sq = ((values - mu) ** 2) / (sigma * sigma)
        log_values = logsumexp(np.log(probs)[None, :] + c[:, None] * scaled_sq[None, :], axis=1)
        log_rhos = invert_alpha_bound_logrho_vector(q_task, alpha_values, log_values)
        vals = log_rhos / math.log(2.0)
        idx = int(np.argmin(vals))
        if float(vals[idx]) < best_log2:
            best_log2 = float(vals[idx])
            best_alpha = float(alpha_values[idx])

    evaluate(alpha_grid)
    for width, count in ((0.5, 81), (0.1, 61)):
        evaluate(np.linspace(max(1.001, best_alpha - width), best_alpha + width, count))
    return best_log2, best_alpha


def prw_gaussian_log2(case: TimingCase, alpha: float, sigma: float) -> float:
    values, normalized, log_base = grouped_prw_weights(case, alpha)
    positive = normalized > 0.0
    values = values[positive]
    normalized = normalized[positive]
    y, dx = integration_grid(values, sigma)
    logpdf = gaussian_logpdf_grid(y, values, sigma)
    log_terms = logsumexp(np.log(normalized)[None, :] + alpha * logpdf, axis=1) / alpha
    return log2(log_base / alpha + float(logsumexp(log_terms)) + math.log(dx))


def alpha_grid() -> np.ndarray:
    return np.unique(
        np.concatenate([
            np.linspace(1.001, 8.0, 80),
            np.linspace(8.25, 64.0, 90),
            np.linspace(68.0, 256.0, 48),
        ])
    )


def required_gaussian_sigma(case: TimingCase, alpha: float, target_log2: float) -> float:
    no_leakage = log2(grouped_prw_weights(case, alpha)[2] / alpha)
    if no_leakage > target_log2:
        return float("inf")
    if prw_gaussian_log2(case, alpha, 0.25) <= target_log2:
        return 0.0
    lo, hi = 0.25, 1.0
    while prw_gaussian_log2(case, alpha, hi) > target_log2:
        hi *= 2.0
        if hi > 1e6:
            return float("inf")
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if prw_gaussian_log2(case, alpha, mid) <= target_log2:
            hi = mid
        else:
            lo = mid
    return hi


def local_prw_target_log2(bits: int, total_target_exponent: float, rounds: int, alpha: float) -> float:
    log_base_bits = (1.0 - alpha) * float(bits)
    total_log_r_bits = alpha * (-float(total_target_exponent)) - log_base_bits
    return (log_base_bits + total_log_r_bits / float(rounds)) / alpha


def map_gaussian_log2(case: TimingCase, sigma: float) -> float:
    values, _ = grouped_probabilities(case)
    y, dx = integration_grid(values, sigma)
    logpdf = gaussian_logpdf_grid(y, values, sigma)
    # Uniform-exponent RSA cases use the same true prior mass for each sampled
    # secret.  The scaled MAP integral is therefore max over timing bins.
    log_scaled = np.max(logpdf, axis=1)
    return case.prior_success_log2 + log2(float(logsumexp(log_scaled)) + math.log(dx))


def required_map_gaussian_sigma(case: TimingCase, target_log2: float) -> float:
    no_leakage = case.prior_success_log2
    if no_leakage > target_log2:
        return float("inf")
    if map_gaussian_log2(case, 0.25) <= target_log2:
        return 0.0
    lo, hi = 0.25, 1.0
    while map_gaussian_log2(case, hi) > target_log2:
        hi *= 2.0
        if hi > 1e6:
            return float("inf")
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if map_gaussian_log2(case, mid) <= target_log2:
            hi = mid
        else:
            lo = mid
    return hi


def scaled_map_value(values: np.ndarray, noise: np.ndarray, probs: np.ndarray) -> float:
    index = {int(v): i for i, v in enumerate(noise)}
    total = 0.0
    for y in range(int(values[0] + noise[0]), int(values[-1] + noise[-1]) + 1):
        best = 0.0
        for h in values:
            j = index.get(int(y - h))
            if j is not None:
                best = max(best, float(probs[j]))
        total += best
    return float(total)


def solve_scaled_map_noise(case: TimingCase, target_log2: float, kind: str, radius: int) -> NoiseResult:
    values = np.unique(case.leakage_bins).astype(int)
    target_scaled = 2.0 ** (target_log2 - case.prior_success_log2)
    if target_scaled < 1.0:
        return NoiseResult("infeasible_below_no_leakage_bound", float("inf"), float("inf"), float("nan"), float("inf"))
    if kind == "symmetric":
        noise = np.arange(-radius, radius + 1, dtype=int)
    elif kind == "onesided":
        noise = np.arange(0, radius + 1, dtype=int)
    else:
        raise ValueError(kind)

    zero = np.zeros(noise.size, dtype=float)
    zero_idx = np.where(noise == 0)[0]
    if zero_idx.size:
        zero[int(zero_idx[0])] = 1.0
        zero_value = scaled_map_value(values, noise, zero)
        if zero_value <= target_scaled + 1e-10:
            return NoiseResult("zero_noise_feasible", 0.0, 0.0, 0.0, zero_value, 0.0)

    n_q = noise.size
    ys = list(range(int(values[0] + noise[0]), int(values[-1] + noise[-1]) + 1))
    n_t = len(ys)
    n_var = n_q + n_t
    index = {int(v): i for i, v in enumerate(noise)}
    rows = []
    rhs = []
    rows.append({n_q + yi: 1.0 for yi in range(n_t)})
    rhs.append(target_scaled)
    for yi, y in enumerate(ys):
        for h in values:
            j = index.get(int(y - h))
            if j is not None:
                rows.append({j: 1.0, n_q + yi: -1.0})
                rhs.append(0.0)
    a_ub = lil_matrix((len(rows), n_var), dtype=float)
    for i, entries in enumerate(rows):
        for j, value in entries.items():
            a_ub[i, j] = value
    a_eq_rows = 1 if kind != "symmetric" else 2
    a_eq = lil_matrix((a_eq_rows, n_var), dtype=float)
    a_eq[0, :n_q] = np.ones(n_q)
    b_eq = [1.0]
    if kind == "symmetric":
        a_eq[1, :n_q] = noise.astype(float)
        b_eq.append(0.0)
    costs = np.abs(noise.astype(float))
    c = np.zeros(n_var, dtype=float)
    c[:n_q] = costs
    opt = linprog(
        c,
        A_ub=a_ub.tocsr(),
        b_ub=np.asarray(rhs, dtype=float),
        A_eq=a_eq.tocsr(),
        b_eq=np.asarray(b_eq, dtype=float),
        bounds=[(0.0, None)] * n_var,
        method="highs",
    )
    if not opt.success or opt.x is None:
        return NoiseResult(str(opt.message), float("inf"), float("inf"), float("nan"), float("inf"))
    probs = np.maximum(opt.x[:n_q], 0.0)
    total = float(np.sum(probs))
    if total > 0:
        probs /= total
    objective = float(np.sum(costs * probs))
    expected_abs = float(np.sum(np.abs(noise.astype(float)) * probs))
    second_moment = float(np.sum((noise.astype(float) ** 2) * probs))
    mean = float(np.sum(noise.astype(float) * probs))
    value = scaled_map_value(values, noise, probs)
    return NoiseResult("optimal", objective, expected_abs, mean, value, second_moment)


def interval_map_noise(case: TimingCase, target_log2: float, kind: str, radius: int) -> NoiseResult:
    """Fast MAP interval optimizer for large-support timing perturbations.

    For additive noise over the enclosing timing-bin interval, the MAP inflation
    factor is 1 + span / support size.  The exact LP hard point at large radius
    matches this value, while a full LP sweep would be prohibitively slow.
    """
    values = np.unique(case.leakage_bins).astype(int)
    span = int(values[-1] - values[0])
    factor = 2.0 ** (target_log2 - case.prior_success_log2)
    if factor <= 1.0:
        return NoiseResult("infeasible_below_no_leakage_bound", float("inf"), float("inf"), float("nan"), float("inf"))
    if factor >= float(values.size):
        return NoiseResult("zero_noise_feasible", 0.0, 0.0, 0.0, float(values.size), 0.0)
    if kind == "symmetric":
        max_support_size = float(2 * radius + 1)
        min_factor = 1.0 + float(span) / max_support_size
        if factor < min_factor:
            return NoiseResult("infeasible_support_radius", float("inf"), float("inf"), float("nan"), float("inf"))
        expected_abs = float(span) / (4.0 * (factor - 1.0))
    elif kind == "onesided":
        max_support_size = float(radius + 1)
        min_factor = 1.0 + float(span) / max_support_size
        if factor < min_factor:
            return NoiseResult("infeasible_support_radius", float("inf"), float("inf"), float("nan"), float("inf"))
        expected_abs = float(span) / (2.0 * (factor - 1.0))
    else:
        raise ValueError(kind)
    return NoiseResult("interval_map_closed_form", expected_abs, expected_abs, float("nan"), factor, float("nan"))


def scaled_prw_value(values: np.ndarray, normalized: np.ndarray, noise: np.ndarray, probs: np.ndarray, alpha: float) -> float:
    index = {int(v): i for i, v in enumerate(noise)}
    total = 0.0
    for y in range(int(values[0] + noise[0]), int(values[-1] + noise[-1]) + 1):
        inner = 0.0
        for h, weight in zip(values, normalized):
            j = index.get(int(y - h))
            if j is not None:
                inner += float(weight) * float(probs[j]) ** alpha
        if inner > 0.0:
            total += inner ** (1.0 / alpha)
    return float(total)


@dataclass
class NoiseResult:
    status: str
    objective: float
    expected_abs: float
    mean: float
    scaled_value: float
    second_moment: float = float("nan")


def solve_scaled_prw_noise(case: TimingCase, alpha: float, target_log2: float, kind: str, radius: int) -> NoiseResult:
    values_f, normalized, log_base = grouped_prw_weights(case, alpha)
    values = values_f.astype(int)
    log_target_scaled = target_log2 * math.log(2.0) - log_base / alpha
    no_leakage_scaled = 1.0
    if math.exp(min(700.0, log_target_scaled)) < no_leakage_scaled:
        return NoiseResult("infeasible_below_no_leakage_bound", float("inf"), float("inf"), float("nan"), float("inf"))
    if kind == "symmetric":
        noise = np.arange(-radius, radius + 1, dtype=int)
    elif kind == "onesided":
        noise = np.arange(0, radius + 1, dtype=int)
    else:
        raise ValueError(kind)
    zero = np.zeros(noise.size, dtype=float)
    zero_idx = np.where(noise == 0)[0]
    if zero_idx.size:
        zero[int(zero_idx[0])] = 1.0
        zero_value = scaled_prw_value(values, normalized, noise, zero, alpha)
        if math.log(zero_value) <= log_target_scaled + 1e-10:
            return NoiseResult("zero_noise_feasible", 0.0, 0.0, 0.0, zero_value, 0.0)

    q = cp.Variable(noise.size, nonneg=True)
    constraints = [cp.sum(q) == 1.0]
    if kind == "symmetric":
        constraints.append(noise.astype(float) @ q == 0.0)
    index = {int(v): i for i, v in enumerate(noise)}
    terms = []
    for y in range(int(values[0] + noise[0]), int(values[-1] + noise[-1]) + 1):
        pieces = []
        for h, weight in zip(values, normalized):
            j = index.get(int(y - h))
            if j is not None and weight > 0.0:
                pieces.append((float(weight) ** (1.0 / alpha)) * q[j])
        if pieces:
            terms.append(cp.norm(cp.hstack(pieces), p=alpha))
    constraints.append(cp.sum(cp.hstack(terms)) <= math.exp(log_target_scaled))
    costs = np.abs(noise.astype(float))
    problem = cp.Problem(cp.Minimize(costs @ q), constraints)
    problem.solve(solver="CLARABEL")
    if q.value is None:
        return NoiseResult(str(problem.status), float("inf"), float("inf"), float("nan"), float("inf"))
    probs = np.maximum(np.asarray(q.value, dtype=float), 0.0)
    total = float(np.sum(probs))
    if total > 0.0:
        probs /= total
    objective = float(np.sum(costs * probs))
    expected_abs = float(np.sum(np.abs(noise.astype(float)) * probs))
    second_moment = float(np.sum((noise.astype(float) ** 2) * probs))
    mean = float(np.sum(noise.astype(float) * probs))
    value = scaled_prw_value(values, normalized, noise, probs, alpha)
    return NoiseResult(str(problem.status), objective, expected_abs, mean, value, second_moment)


def curves_for_case(case: TimingCase, sigmas: np.ndarray, rows: List[dict]) -> Dict[str, List[float]]:
    pac_grid = alpha_grid()
    curves = {
        "Fano Mutual Information": [],
        "PAC-alpha(best)": [],
        "PRW alpha=16": [],
        "PRW alpha=64": [],
        "PRW alpha=256": [],
    }
    q_task = 2.0 ** case.prior_success_log2
    for sigma in sigmas:
        mi = mutual_information(case, float(sigma))
        pac_kl = math.log2(invert_kl_bound(q_task, mi))
        pac_alpha, best_alpha = pac_alpha_best_log2(case, float(sigma), pac_grid)
        vals = [pac_kl, pac_alpha]
        for alpha in PRW_ALPHAS:
            vals.append(min(0.0, prw_gaussian_log2(case, alpha, float(sigma)) + case.task_log2_bonus))
        for key, value in zip(curves.keys(), vals):
            curves[key].append(value)
        rows.append({
            "panel": case.label,
            "sigma": float(sigma),
            "fano_mi_log2": pac_kl,
            "pac_alpha_best_log2": pac_alpha,
            "pac_alpha_best_alpha": best_alpha,
            "prw_alpha_16_log2": vals[2],
            "prw_alpha_64_log2": vals[3],
            "prw_alpha_256_log2": vals[4],
            "task_log2_bonus": case.task_log2_bonus,
            "prior_success_log2": case.prior_success_log2,
        })
    return curves


def build_results() -> Tuple[List[dict], List[dict]]:
    rows: List[dict] = []
    panels: List[dict] = []
    sigmas = np.arange(0.5, 80.0 + 1e-9, 1.0)

    case128 = make_timing_case(128, label="(a) Uniform 128-bit, full identification", seed=SEED + 1)
    case256 = make_timing_case(256, label="(b) Uniform 256-bit, full identification", seed=SEED + 2)
    case256_poisson = make_poisson_hw_timing_case(
        256,
        label="(c) Poisson-HW 256-bit full identification",
        seed=SEED + 3,
        lam=96.0,
    )
    case256_recover240 = make_timing_case(
        256,
        label="(d) Uniform 256-bit, recover at least 240 bit",
        seed=SEED + 4,
    )
    recover_bonus = log2_hamming_ball_size(256, 16)
    case256_recover240.task_log2_bonus = recover_bonus
    case256_recover240.prior_success_log2 = -256.0 + recover_bonus
    base256 = case256

    for case in [case128, case256, case256_poisson, case256_recover240]:
        panels.append({
            "title": case.label,
            "series": curves_for_case(case, sigmas, rows),
            "x_values": list(map(float, sigmas)),
            "x_ticks": [0.5, 20, 40, 60, 80],
            "ylabel": "",
            "kind": "comparison",
        })

    target_exponents_e = np.arange(248.0, 252.0 + 1e-9, 1.0)
    entropy_series: Dict[str, List[Tuple[float, float]]] = {}
    for bits in [253, 254, 255, 256]:
        case = make_timing_case(bits, label=f"l={bits}", seed=SEED + 10 + bits)
        pts = []
        for exponent in target_exponents_e:
            sigma = required_map_gaussian_sigma(case, -float(exponent))
            pts.append((float(exponent), sigma))
            rows.append({
                "panel": "(e) entropy",
                "secret_bits": bits,
                "target_exponent": exponent,
                "sigma": sigma,
                "alpha": OPT_ALPHA,
            })
        entropy_series[f"l={bits}"] = pts
    panels.append({
        "title": "(e) Effect from various $l$-bit secret entropy",
        "series": entropy_series,
        "x_ticks": [248, 249, 250, 251, 252],
        "ylabel": "Gaussian noise std σ",
        "kind": "entropy",
    })

    target_exponents_f = np.arange(250.0, 254.0 + 1e-9, 1.0)
    mechanism_series: Dict[str, List[Tuple[float, float]]] = {
        "Gaussian": [],
        BOUNDED_LABEL: [],
        ONESIDED_LABEL: [],
    }
    for exponent in target_exponents_f:
        target = -float(exponent)
        sigma = required_map_gaussian_sigma(base256, target)
        mechanism_series["Gaussian"].append((float(exponent), sigma * math.sqrt(2.0 / math.pi)))
        for label, kind in [
            (BOUNDED_LABEL, "symmetric"),
            (ONESIDED_LABEL, "onesided"),
        ]:
            result = solve_scaled_map_noise(base256, target, kind, NOISE_RADIUS)
            mechanism_series[label].append((float(exponent), result.expected_abs))
            rows.append({
                "panel": "(f) mechanisms",
                "target_exponent": exponent,
                "mechanism": label,
                "expected_abs": result.expected_abs,
                "objective_value": result.objective,
                "optimization_objective": "expected_abs",
                "second_moment": result.second_moment,
                "mean": result.mean,
                "scaled_prw_value": result.scaled_value,
                "status": result.status,
                "noise_radius": NOISE_RADIUS,
                "accounting": "map_alpha_infinity_lp",
            })
    panels.append({
        "title": "(f) Randomization mechanisms",
        "series": mechanism_series,
        "x_ticks": [250, 251, 252, 253, 254],
        "ylabel": "Expected absolute perturbation",
        "kind": "mechanism",
    })

    comp_rng = np.random.default_rng(SEED + 1000)
    composition_trials: List[TimingCase] = []
    for trial in range(COMPOSITION_TRIALS):
        query = make_fast_composition_query(comp_rng, base256.bits)
        composition_trials.append(
            make_fast_composition_case(
                base256.samples,
                bits=base256.bits,
                label=f"composition representative trial {trial + 1}",
                query=query,
                log_prior=base256.log_prior,
                prior_success_log2=base256.prior_success_log2,
            )
        )

    for letter, total_exponent in [("g", 240.0), ("h", 220.0)]:
        comp_series: Dict[str, List[Tuple[float, float]]] = {
            "Gaussian": [],
            BOUNDED_LABEL: [],
            ONESIDED_LABEL: [],
        }
        for rounds in COMPOSITION_ROUNDS:
            local_log2_delta = base256.prior_success_log2 + (base256.bits - total_exponent) / rounds
            local_prw_log2 = local_prw_target_log2(
                base256.bits,
                total_exponent,
                rounds,
                COMPOSITION_PRW_ALPHA,
            )

            gaussian_trial_values = []
            for trial_idx, case in enumerate(composition_trials, start=1):
                sigma = required_gaussian_sigma(case, COMPOSITION_PRW_ALPHA, local_prw_log2)
                gaussian_abs = sigma * math.sqrt(2.0 / math.pi)
                gaussian_trial_values.append(float(gaussian_abs))
                rows.append({
                    "panel": f"({letter}) composition",
                    "target_exponent": total_exponent,
                    "rounds": rounds,
                    "trial": trial_idx,
                    "mechanism": "Gaussian",
                    "local_log2_delta": local_log2_delta,
                    "local_prw_target_log2": local_prw_log2,
                    "expected_abs": gaussian_abs,
                    "expected_abs_std_within_trial": 0.0,
                    "reported_quantity": "prefix_average_expected_abs_noise",
                    "status": "finite_alpha_prw_equal_split",
                    "composition_alpha": COMPOSITION_PRW_ALPHA,
                })

            comp_series["Gaussian"].append((float(rounds), float(np.mean(gaussian_trial_values))))
            rows.append({
                "panel": f"({letter}) composition",
                "target_exponent": total_exponent,
                "rounds": rounds,
                "mechanism": "Gaussian",
                "local_log2_delta": local_log2_delta,
                "local_prw_target_log2": local_prw_log2,
                "expected_abs": float(np.mean(gaussian_trial_values)),
                "expected_abs_across_trials_std": float(np.std(gaussian_trial_values)),
                "trials": COMPOSITION_TRIALS,
                "reported_quantity": "prefix_average_expected_abs_noise",
                "status": "finite_alpha_prw_equal_split",
                "composition_alpha": COMPOSITION_PRW_ALPHA,
            })

            for label, kind in [
                (BOUNDED_LABEL, "symmetric"),
                (ONESIDED_LABEL, "onesided"),
            ]:
                trial_values = []
                trial_statuses = []
                for trial_idx, case in enumerate(composition_trials, start=1):
                    result = solve_scaled_map_noise(case, local_log2_delta, kind, NOISE_RADIUS)
                    trial_values.append(float(result.expected_abs))
                    trial_statuses.append(result.status)
                    rows.append({
                        "panel": f"({letter}) composition",
                        "target_exponent": total_exponent,
                        "rounds": rounds,
                        "trial": trial_idx,
                        "mechanism": label,
                        "local_log2_delta": local_log2_delta,
                        "local_prw_target_log2": local_prw_log2,
                        "expected_abs": result.expected_abs,
                        "reported_quantity": "prefix_average_expected_abs_noise",
                        "status": result.status,
                        "noise_radius": NOISE_RADIUS,
                        "accounting": "map_alpha_infinity_lp",
                    })
                finite_values = [v for v in trial_values if math.isfinite(v)]
                mean_value = float(np.mean(finite_values)) if len(finite_values) == len(trial_values) else float("inf")
                comp_series[label].append((float(rounds), mean_value))
                rows.append({
                    "panel": f"({letter}) composition",
                    "target_exponent": total_exponent,
                    "rounds": rounds,
                    "mechanism": label,
                    "local_log2_delta": local_log2_delta,
                    "local_prw_target_log2": local_prw_log2,
                    "expected_abs": mean_value,
                    "expected_abs_across_trials_std": float(np.std(finite_values)) if len(finite_values) == len(trial_values) else "",
                    "trials": COMPOSITION_TRIALS,
                    "reported_quantity": "prefix_average_expected_abs_noise",
                    "status": "all_trials_feasible" if len(finite_values) == len(trial_values) else "some_trials_infeasible",
                    "trial_statuses": ";".join(trial_statuses),
                    "noise_radius": NOISE_RADIUS,
                    "accounting": "map_alpha_infinity_lp",
                })
        panels.append({
            "title": f"({letter}) Composition, full identification $2^{{-{int(total_exponent)}}}$",
            "series": comp_series,
            "x_ticks": COMPOSITION_X_TICKS,
            "ylabel": "",
            "kind": "composition",
        })

    meta = [
        {"key": "sample_count", "value": SAMPLE_COUNT},
        {"key": "seed", "value": SEED},
        {"key": "timing_bin_unit", "value": "1 ns timing bins"},
        {"key": "ground_truth", "value": "not plotted; RSA timing is empirical black-box"},
        {"key": "optimization_alpha", "value": OPT_ALPHA},
        {"key": "noise_radius", "value": NOISE_RADIUS},
        {"key": "composition_trials", "value": COMPOSITION_TRIALS},
        {"key": "composition_query_model", "value": "representative random vectorized RSA-style timing query per trial"},
        {"key": "composition_reported_quantity", "value": "prefix average t^{-1} sum_{i=1}^t E|N_i|"},
        {"key": "composition_alpha", "value": COMPOSITION_PRW_ALPHA},
    ]
    return panels, rows + meta


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    fields: List[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def draw_rich_text(d: ImageDraw.ImageDraw, xy: Tuple[float, float], text: str, text_font, *, scale: int) -> None:
    x, y = xy
    small = font(max(1, int(text_font.size * 0.78)), True)
    i = 0
    while i < len(text):
        if text.startswith("$l$", i):
            d.text((x, y), "l", fill=(25, 25, 25), font=text_font)
            x += d.textlength("l", font=text_font)
            i += 3
            continue
        if text.startswith("$2^{-", i):
            end = text.find("}$", i)
            if end != -1:
                exponent = "-" + text[i + 5:end]
                d.text((x, y), "2", fill=(25, 25, 25), font=text_font)
                x += d.textlength("2", font=text_font) + scale
                d.text((x, y - 12 * scale), exponent, fill=(25, 25, 25), font=small)
                x += d.textlength(exponent, font=small)
                i = end + 2
                continue
        next_math = text.find("$", i)
        chunk = text[i:] if next_math == -1 else text[i:next_math]
        d.text((x, y), chunk, fill=(25, 25, 25), font=text_font)
        x += d.textlength(chunk, font=text_font)
        i += len(chunk)


def draw_log2_label(d: ImageDraw.ImageDraw, center_x: float, y: float, suffix: str, axis_font, *, scale: int) -> None:
    sub = font(max(1, int(axis_font.size * 0.78)), True)
    total = d.textlength("log", font=axis_font) + d.textlength("2", font=sub) + d.textlength(suffix, font=axis_font) + 11 * scale
    x = center_x - total / 2
    d.text((x, y), "log", fill=(25, 25, 25), font=axis_font)
    x += d.textlength("log", font=axis_font) + 3 * scale
    d.text((x, y + 26 * scale), "2", fill=(25, 25, 25), font=sub)
    x += d.textlength("2", font=sub) + 8 * scale
    d.text((x, y), suffix, fill=(25, 25, 25), font=axis_font)


def draw_rotated_label(img: Image.Image, text: str, center_x: float, center_y: float, text_font, *, scale: int) -> None:
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    w = int(measure.textlength(text, font=text_font) + 24 * scale)
    tmp = Image.new("RGBA", (w, 70 * scale), (255, 255, 255, 0))
    td = ImageDraw.Draw(tmp)
    td.text((0, 0), text, fill=(25, 25, 25), font=text_font)
    bbox = tmp.getbbox()
    if bbox:
        tmp = tmp.crop(bbox)
    rot = tmp.rotate(90, expand=True)
    img.paste(rot, (int(center_x - rot.width / 2), int(center_y - rot.height / 2)), rot)


def color_map(name: str, kind: str) -> Tuple[int, int, int]:
    comparison = {
        "Fano Mutual Information": (0, 114, 178),
        "PAC-alpha(best)": (213, 94, 0),
        "PRW alpha=16": (0, 158, 115),
        "PRW alpha=64": (230, 159, 0),
        "PRW alpha=256": (204, 121, 167),
        "PRW required sigma": (0, 158, 115),
    }
    entropy = {
        "l=253": (0, 114, 178),
        "l=254": (213, 94, 0),
        "l=255": (0, 158, 115),
        "l=256": (204, 121, 167),
    }
    mechanism = {
        "Gaussian": (86, 180, 233),
        BOUNDED_LABEL: (117, 112, 179),
        ONESIDED_LABEL: (166, 86, 40),
    }
    if kind == "entropy":
        return entropy[name]
    if kind in {"mechanism", "composition"}:
        return mechanism[name]
    return comparison[name]


def y_label(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2g}"


def draw_comparison_zoom_inset(
    d: ImageDraw.ImageDraw,
    series: Dict[str, List[Tuple[float, float]]],
    panel_title: str,
    main_polylines: Sequence[Sequence[Tuple[float, float]]],
    box: Tuple[float, float, float, float],
    tick_font,
    *,
    scale: int,
) -> None:
    zoom_names = [
        "PAC-alpha(best)",
        "PRW alpha=16",
        "PRW alpha=64",
        "PRW alpha=256",
    ]
    zoom_series = {
        name: [(x, y) for x, y in series.get(name, []) if 40.0 <= x <= 80.0 and math.isfinite(y)]
        for name in zoom_names
    }
    zoom_series = {name: pts for name, pts in zoom_series.items() if pts}
    if not zoom_series:
        return

    x0, y0, x1, y1 = box
    w = x1 - x0
    h = y1 - y0

    def orient(a: Tuple[float, float], b: Tuple[float, float], c: Tuple[float, float]) -> float:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def segments_intersect(
        a: Tuple[float, float],
        b: Tuple[float, float],
        c: Tuple[float, float],
        e: Tuple[float, float],
    ) -> bool:
        o1 = orient(a, b, c)
        o2 = orient(a, b, e)
        o3 = orient(c, e, a)
        o4 = orient(c, e, b)
        return (o1 == 0.0 or o2 == 0.0 or o1 * o2 < 0.0) and (o3 == 0.0 or o4 == 0.0 or o3 * o4 < 0.0)

    def segment_hits_rect(a: Tuple[float, float], b: Tuple[float, float], rect: Tuple[float, float, float, float]) -> bool:
        left, top, right, bottom = rect
        if max(a[0], b[0]) < left or min(a[0], b[0]) > right or max(a[1], b[1]) < top or min(a[1], b[1]) > bottom:
            return False
        if left <= a[0] <= right and top <= a[1] <= bottom:
            return True
        if left <= b[0] <= right and top <= b[1] <= bottom:
            return True
        corners = [(left, top), (right, top), (right, bottom), (left, bottom)]
        edges = list(zip(corners, corners[1:] + corners[:1]))
        return any(segments_intersect(a, b, edge_a, edge_b) for edge_a, edge_b in edges)

    def rect_hits_curves(rect: Tuple[float, float, float, float], pad: float) -> bool:
        padded = (rect[0] - pad, rect[1] - pad, rect[2] + pad, rect[3] + pad)
        for pts in main_polylines:
            for a, b in zip(pts, pts[1:]):
                if segment_hits_rect(a, b, padded):
                    return True
        return False

    def choose_rect() -> Tuple[float, float, float, float]:
        preferred_y = 0.27 if panel_title.startswith("(c)") else 0.48
        preferred = (x0 + 0.70 * w, y0 + preferred_y * h)
        best: Tuple[float, float, float, float] | None = None
        best_score = float("inf")
        base_h = 0.34 if panel_title.startswith("(c)") else 0.42
        for factor in [1.0, 0.94, 0.88, 0.82, 0.76, 0.70]:
            iw = 0.60 * w * factor
            ih = base_h * h * factor
            x_candidates = np.linspace(x0 + 0.04 * w, x1 - iw - 0.02 * w, 10)
            y_candidates = np.linspace(y0 + 0.08 * h, y1 - ih - 0.05 * h, 14)
            for ix0_candidate in x_candidates:
                for iy0_candidate in y_candidates:
                    rect = (float(ix0_candidate), float(iy0_candidate), float(ix0_candidate + iw), float(iy0_candidate + ih))
                    if rect_hits_curves(rect, 9 * scale):
                        continue
                    center = ((rect[0] + rect[2]) / 2.0, (rect[1] + rect[3]) / 2.0)
                    score = abs(center[0] - preferred[0]) + 0.75 * abs(center[1] - preferred[1])
                    if score < best_score:
                        best_score = score
                        best = rect
            if best is not None:
                return best
        return (x1 - 0.42 * w, y0 + 0.30 * h, x1 - 0.02 * w, y0 + 0.56 * h)

    ix0, iy0, ix1, iy1 = choose_rect()
    iw = ix1 - ix0
    ih = iy1 - iy0
    xmin, xmax = 40.0, 80.0
    all_y = [y for pts in zoom_series.values() for _, y in pts]
    ymin = min(all_y)
    ymax = max(all_y)
    pad = max(0.4, 0.08 * (ymax - ymin))
    ymin = math.floor(ymin - pad)
    ymax = math.ceil(ymax + pad)
    if ymax <= ymin:
        ymax = ymin + 1.0

    def xm(value: float) -> float:
        return ix0 + (value - xmin) / (xmax - xmin) * iw

    def ym(value: float) -> float:
        return iy0 + (ymax - value) / (ymax - ymin) * ih

    inset_title_font = font(25 * scale, True)
    inset_tick_font = font(24 * scale)
    d.rectangle([ix0, iy0, ix1, iy1], fill=(255, 255, 255), outline=(120, 120, 120), width=2 * scale)
    for tick in [40, 60, 80]:
        px = xm(float(tick))
        d.line([(px, iy0), (px, iy1)], fill=(232, 232, 232), width=1 * scale)
        label = f"{tick:g}"
        tw = d.textlength(label, font=inset_tick_font)
        d.text((px - tw / 2, iy1 - 31 * scale), label, fill=(55, 55, 55), font=inset_tick_font)
    for tick in np.linspace(ymin, ymax, 3):
        py = ym(float(tick))
        d.line([(ix0, py), (ix1, py)], fill=(232, 232, 232), width=1 * scale)
        label = y_label(float(tick))
        d.text((ix0 + 8 * scale, py - 13 * scale), label, fill=(55, 55, 55), font=inset_tick_font)
    d.text((ix0 + 10 * scale, iy0 + 7 * scale), "σ=40-80", fill=(35, 35, 35), font=inset_title_font)

    for name, pts in zoom_series.items():
        color = color_map(name, "comparison")
        xy = [(xm(x), ym(y)) for x, y in pts]
        for a, b in zip(xy, xy[1:]):
            d.line([a, b], fill=color, width=3 * scale)


def draw_panel(
    img: Image.Image,
    d: ImageDraw.ImageDraw,
    panel: dict,
    box: Tuple[float, float, float, float],
    fonts: Tuple[object, object, object],
    *,
    scale: int,
) -> None:
    title_font, tick_font, legend_font = fonts
    x0, y0, x1, y1 = box
    w = x1 - x0
    h = y1 - y0
    kind = str(panel["kind"])
    raw_series = panel["series"]
    if kind == "comparison":
        raw_series = {name: pts for name, pts in raw_series.items() if name != "Ground truth"}
    x_values = panel.get("x_values")
    series: Dict[str, List[Tuple[float, float]]] = {}
    for name, pts in raw_series.items():
        if pts and isinstance(pts[0], tuple):
            series[name] = pts
        elif x_values is not None:
            series[name] = [(float(x), float(y)) for x, y in zip(x_values, pts)]
        else:
            raise TypeError(f"series {name} needs point pairs or panel x_values")
    all_x = [x for pts in series.values() for x, y in pts if math.isfinite(y)]
    all_y = [y for pts in series.values() for x, y in pts if math.isfinite(y)]
    xmin, xmax = min(all_x), max(all_x)
    ymin = min(0.0, min(all_y)) if kind in {"entropy", "mechanism", "composition", "alpha"} else math.floor(min(all_y) - 1.0)
    ymax = max(all_y)
    ymax = ymin + 1.12 * (ymax - ymin) if ymin >= 0.0 else math.ceil(ymax + 1.0)
    if ymax <= ymin:
        ymax = ymin + 1.0

    def xm(value: float) -> float:
        return x0 + (value - xmin) / (xmax - xmin) * w

    def ym(value: float) -> float:
        return y0 + (ymax - value) / (ymax - ymin) * h

    draw_rich_text(d, (x0, y0 - 58 * scale), str(panel["title"]), title_font, scale=scale)
    d.rectangle([x0, y0, x1, y1], fill=(250, 250, 250), outline=(178, 178, 178))

    for tick in panel["x_ticks"]:
        px = xm(float(tick))
        d.line([(px, y0), (px, y1)], fill=(228, 228, 228))
        label = f"{tick:g}"
        tw = d.textlength(label, font=tick_font)
        d.text((px - tw / 2, y1 + 12 * scale), label, fill=(45, 45, 45), font=tick_font)

    for idx, tick in enumerate(np.linspace(ymin, ymax, 6)):
        py = ym(float(tick))
        d.line([(x0, py), (x1, py)], fill=(228, 228, 228))
        if idx == 0:
            continue
        label = y_label(float(tick))
        tw = d.textlength(label, font=tick_font)
        d.text((x0 - tw - 10 * scale, py - 13 * scale), label, fill=(45, 45, 45), font=tick_font)

    d.line([(x0, y1), (x1, y1)], fill=(30, 30, 30), width=3 * scale)
    d.line([(x0, y0), (x0, y1)], fill=(30, 30, 30), width=3 * scale)

    main_polylines: List[List[Tuple[float, float]]] = []
    for name, pts in series.items():
        color = color_map(name, kind)
        clean = [(x, y) for x, y in pts if math.isfinite(y)]
        xy = [(xm(x), ym(y)) for x, y in clean]
        main_polylines.append(xy)
        for a, b in zip(xy, xy[1:]):
            d.line([a, b], fill=color, width=5 * scale)
        if kind in {"alpha", "composition"}:
            for px, py in xy:
                d.ellipse([px - 4 * scale, py - 4 * scale, px + 4 * scale, py + 4 * scale], fill=color, outline="white")

    if kind == "comparison":
        draw_comparison_zoom_inset(d, series, str(panel["title"]), main_polylines, box, tick_font, scale=scale)

    if kind == "entropy":
        lx, ly = x0 + 22 * scale, y0 + 20 * scale
        for idx, name in enumerate(series.keys()):
            yy = ly + idx * 64 * scale
            color = color_map(name, kind)
            d.line([(lx, yy + 16 * scale), (lx + 42 * scale, yy + 16 * scale)], fill=color, width=5 * scale)
            d.text((lx + 52 * scale, yy), name, fill=(28, 28, 28), font=legend_font)


def draw_legend(
    d: ImageDraw.ImageDraw,
    labels: Sequence[str],
    kind: str,
    center_x: float,
    y: float,
    legend_font,
    *,
    scale: int,
) -> None:
    widths = [52 * scale + d.textlength(label.replace("alpha", "α"), font=legend_font) + 30 * scale for label in labels]
    x = center_x - sum(widths) / 2.0
    for label in labels:
        text = label.replace("alpha", "α")
        color = color_map(label, kind)
        d.line([(x, y + 16 * scale), (x + 42 * scale, y + 16 * scale)], fill=color, width=5 * scale)
        d.text((x + 52 * scale, y), text, fill=(28, 28, 28), font=legend_font)
        x += 52 * scale + d.textlength(text, font=legend_font) + 30 * scale


def plot_combined(path: Path, panels: Sequence[dict]) -> None:
    final_width, final_height = 5600, 2500
    scale = 2
    width, height = final_width * scale, final_height * scale
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    title_font = font(56 * scale, True)
    tick_font = font(56 * scale)
    axis_font = font(56 * scale, True)
    legend_font = font(56 * scale)
    left, right, gap = 260 * scale, 60 * scale, 220 * scale
    panel_w = (width - left - right - 3 * gap) / 4.0
    panel_h = 850 * scale
    top_y, bottom_y = 190 * scale, 1450 * scale

    comparison_labels = [name for name in panels[0]["series"].keys() if name != "Ground truth"]
    draw_legend(d, comparison_labels, "comparison", width / 2.0, 38 * scale, legend_font, scale=scale)
    for row, y0 in enumerate([top_y, bottom_y]):
        for col in range(4):
            idx = row * 4 + col
            x0 = left + col * (panel_w + gap)
            draw_panel(
                img,
                d,
                panels[idx],
                (x0, y0, x0 + panel_w, y0 + panel_h),
                (title_font, tick_font, legend_font),
                scale=scale,
            )

    draw_rotated_label(img, "log2 of posterior success rate", 60 * scale, top_y + panel_h / 2, axis_font, scale=scale)
    tw = d.textlength("Gaussian timing-noise std σ (timing bins)", font=axis_font)
    d.text((width / 2 - tw / 2, 1140 * scale), "Gaussian timing-noise std σ (timing bins)", fill=(25, 25, 25), font=axis_font)

    f_x0 = left + 1 * (panel_w + gap)
    h_x1 = left + 3 * (panel_w + gap) + panel_w
    draw_legend(d, list(panels[5]["series"].keys()), "mechanism", (f_x0 + h_x1) / 2.0, 1220 * scale, legend_font, scale=scale)

    ef_center = left + (2 * panel_w + gap) / 2.0
    gh_left = left + 2 * (panel_w + gap)
    gh_center = gh_left + (2 * panel_w + gap) / 2.0
    draw_log2_label(d, ef_center, 2405 * scale, " of target posterior success rate", axis_font, scale=scale)
    comp = "Composition iteration"
    tw = d.textlength(comp, font=axis_font)
    d.text((gh_center - tw / 2, 2405 * scale), comp, fill=(25, 25, 25), font=axis_font)
    draw_rotated_label(img, "Gaussian noise std σ", left - 200 * scale, bottom_y + panel_h / 2, axis_font, scale=scale)
    draw_rotated_label(img, "Expected absolute perturbation", f_x0 - 160 * scale, bottom_y + panel_h / 2, axis_font, scale=scale)

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    img = img.resize((final_width, final_height), resampling)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panels, rows = build_results()
    write_csv(OUT / "rsa_timing_empirical_results.csv", rows)
    plot_combined(OUT / "rsa_timing_empirical_combined_v2.png", panels)
    plot_combined(OUT / "rsa_timing_experiments_2x4.png", panels)
    plot_combined(OUT / "rsa_timing_empirical_combined_v8_no_line_overlap_magnifiers.png", panels)
    print(f"wrote {OUT / 'rsa_timing_empirical_results.csv'}")
    print(f"wrote {OUT / 'rsa_timing_empirical_combined_v2.png'}")
    print(f"wrote {OUT / 'rsa_timing_experiments_2x4.png'}")
    print(f"wrote {OUT / 'rsa_timing_empirical_combined_v8_no_line_overlap_magnifiers.png'}")


if __name__ == "__main__":
    main()
