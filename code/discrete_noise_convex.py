"""Convex finite-sequence optimizers for discrete additive noise.

This module implements the finite-domain version of the noise-calibration
problem directly.  Once the additive noise distribution is represented by a
probability vector q over a finite integer support, exact MAP calibration is a
linear program, while finite-alpha PRW calibration is a convex norm-constrained
program.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence

import cvxpy as cp
import numpy as np
from scipy.optimize import linprog
from scipy.special import logsumexp
from scipy.sparse import lil_matrix

from discrete_leakage_experiments import (
    Channel,
    make_aes_power_channel,
    make_rsa_timing_bit_channel,
    n_choose_k_counts,
)


@dataclass
class SequenceNoiseResult:
    status: str
    method: str
    objective: float
    mean: float
    support: np.ndarray
    probabilities: np.ndarray
    risk_per_coord: float
    target_per_coord: float
    full_log_risk: float
    solver: str


def utility_costs(noise: np.ndarray, objective: str) -> np.ndarray:
    if objective == "second_moment":
        return noise.astype(float) ** 2
    if objective == "expected_abs":
        return np.abs(noise.astype(float))
    raise ValueError(f"unknown objective: {objective}")


def _log_sum_prior_alpha(channel: Channel, alpha: float) -> np.ndarray:
    provider = getattr(channel.sum_prior_alpha, "log_values", None)
    if provider is not None:
        return np.asarray(provider(alpha), dtype=float)
    s_alpha = np.asarray(channel.sum_prior_alpha(alpha), dtype=float)
    out = np.full(s_alpha.shape, -math.inf, dtype=float)
    positive = s_alpha > 0.0
    out[positive] = np.log(s_alpha[positive])
    return out


def make_empirical_channel(
    leakage_samples: Sequence[int],
    weights: Optional[Sequence[float]] = None,
    *,
    name: str = "empirical channel",
    coord_bits: int = 1,
    unit_label: str = "bin",
) -> Channel:
    """Build a finite channel from sampled, already-discretized leakages.

    The empirical distribution is supported on the sampled atoms.  If weights
    are omitted, each sampled atom has mass 1/m.  Samples with the same leakage
    value are grouped only for computation; the posterior identification task
    still treats the sampled atoms as the secrets being reconstructed.
    """
    samples = np.asarray(leakage_samples, dtype=int).reshape(-1)
    if samples.size == 0:
        raise ValueError("at least one leakage sample is required")

    if weights is None:
        atom_weights = np.full(samples.size, 1.0 / samples.size, dtype=float)
    else:
        atom_weights = np.asarray(weights, dtype=float).reshape(-1)
        if atom_weights.shape != samples.shape:
            raise ValueError("weights must have the same length as leakage_samples")
        if np.any(atom_weights < 0.0):
            raise ValueError("weights must be nonnegative")
        total = float(np.sum(atom_weights))
        if total <= 0.0:
            raise ValueError("weights must have positive total mass")
        atom_weights = atom_weights / total

    values, inverse = np.unique(samples, return_inverse=True)
    leakage_probs = np.zeros(values.size, dtype=float)
    max_prior_by_leakage = np.zeros(values.size, dtype=float)
    grouped_weights = []
    for i in range(values.size):
        group = atom_weights[inverse == i]
        leakage_probs[i] = float(np.sum(group))
        max_prior_by_leakage[i] = float(np.max(group))
        grouped_weights.append(group.copy())

    def sum_prior_alpha(alpha: float) -> np.ndarray:
        return np.array([np.sum(group**alpha) for group in grouped_weights], dtype=float)

    return Channel(
        name=name,
        leakage_values=values.astype(int),
        leakage_probs=leakage_probs,
        sum_prior_alpha=sum_prior_alpha,
        max_prior_by_leakage=max_prior_by_leakage,
        prior_success_coord=float(np.max(atom_weights)),
        coord_bits=coord_bits,
        unit_label=unit_label,
    )


def make_binomial_aggregate_channel(
    bits: int,
    *,
    bit_prob_one: float = 0.5,
    name: Optional[str] = None,
    unit_label: str = "power bin",
) -> Channel:
    """Analytic aggregate Hamming-weight channel H ~ Bin(bits, bit_prob_one).

    This is equivalent to the complete empirical distribution over all bit
    strings under an independent Bernoulli prior, grouped by their
    one-dimensional Hamming-weight leakage.
    """
    if bits <= 0:
        raise ValueError("bits must be positive")
    p = float(bit_prob_one)
    if not 0.0 < p < 1.0:
        raise ValueError("bit_prob_one must be strictly between zero and one")
    values = np.arange(bits + 1, dtype=int)
    counts = n_choose_k_counts(bits)
    log_atom = values * math.log(p) + (bits - values) * math.log1p(-p)
    probs = counts * np.exp(log_atom)
    probs /= float(np.sum(probs))

    def sum_prior_alpha(alpha: float) -> np.ndarray:
        logs = log_sum_prior_alpha(alpha)
        return np.exp(logs)

    def log_sum_prior_alpha(alpha: float) -> np.ndarray:
        return np.log(counts) + alpha * log_atom

    sum_prior_alpha.log_values = log_sum_prior_alpha  # type: ignore[attr-defined]

    return Channel(
        name=name or f"aggregate HW({bits}, p={p:g})",
        leakage_values=values,
        leakage_probs=probs,
        sum_prior_alpha=sum_prior_alpha,
        max_prior_by_leakage=np.exp(log_atom),
        prior_success_coord=max(p, 1.0 - p) ** bits,
        coord_bits=bits,
        unit_label=unit_label,
    )


def noise_support(kind: str, radius: int) -> np.ndarray:
    if kind == "symmetric":
        return np.arange(-radius, radius + 1, dtype=int)
    if kind == "onesided":
        return np.arange(0, radius + 1, dtype=int)
    if kind in {"nonpositive", "negative_onesided"}:
        return np.arange(-radius, 1, dtype=int)
    raise ValueError(f"unknown support kind: {kind}")


def output_support(channel: Channel, noise: np.ndarray) -> range:
    lo = int(channel.leakage_values[0] + noise[0])
    hi = int(channel.leakage_values[-1] + noise[-1])
    return range(lo, hi + 1)


def _clean_probabilities(q_value: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    q = np.asarray(q_value, dtype=float).reshape(-1)
    q[np.abs(q) < tol] = 0.0
    q = np.maximum(q, 0.0)
    total = float(np.sum(q))
    if total > 0.0:
        q /= total
    return q


def _result(
    status: str,
    method: str,
    objective: float,
    support: np.ndarray,
    probabilities: np.ndarray,
    risk_per_coord: float,
    target_per_coord: float,
    coord_count: int,
    solver: str,
) -> SequenceNoiseResult:
    mean = float(np.sum(support.astype(float) * probabilities))
    full_log_risk = coord_count * math.log(max(risk_per_coord, 1e-300))
    return SequenceNoiseResult(
        status=status,
        method=method,
        objective=float(objective),
        mean=mean,
        support=support,
        probabilities=probabilities,
        risk_per_coord=float(risk_per_coord),
        target_per_coord=float(target_per_coord),
        full_log_risk=float(full_log_risk),
        solver=solver,
    )


def exact_map_risk(channel: Channel, noise: np.ndarray, q: np.ndarray) -> float:
    """Exact one-coordinate MAP identification success."""
    index = {int(v): i for i, v in enumerate(noise)}
    total = 0.0
    for y in output_support(channel, noise):
        best = 0.0
        for h, max_prior in zip(channel.leakage_values, channel.max_prior_by_leakage):
            j = index.get(int(y - h))
            if j is not None:
                best = max(best, float(max_prior) * float(q[j]))
        total += best
    return float(total)


def prw_alpha_value(channel: Channel, noise: np.ndarray, q: np.ndarray, alpha: float) -> float:
    """One-coordinate PRW alpha value on the finite output domain."""
    log_value = log_prw_alpha_value(channel, noise, q, alpha)
    if log_value <= -745.0:
        return 0.0
    return float(math.exp(log_value))


def log_prw_alpha_value(channel: Channel, noise: np.ndarray, q: np.ndarray, alpha: float) -> float:
    """Log one-coordinate PRW alpha value, evaluated stably for large alpha."""
    index = {int(v): i for i, v in enumerate(noise)}
    log_s_alpha = _log_sum_prior_alpha(channel, alpha)
    log_q = np.full(q.shape, -math.inf, dtype=float)
    positive_q = q > 0.0
    log_q[positive_q] = np.log(q[positive_q])
    terms = []
    for y in output_support(channel, noise):
        pieces = []
        for h, log_weight in zip(channel.leakage_values, log_s_alpha):
            j = index.get(int(y - h))
            if j is not None and math.isfinite(float(log_weight)) and math.isfinite(float(log_q[j])):
                pieces.append(float(log_weight) + alpha * float(log_q[j]))
        if pieces:
            terms.append(float(logsumexp(pieces)) / alpha)
    if not terms:
        return -math.inf
    return float(logsumexp(terms))


def solve_exact_map_noise(
    channel: Channel,
    coord_count: int,
    target_log2: float,
    kind: str,
    radius: int,
    *,
    solver: str = "HIGHS",
    objective: str = "second_moment",
) -> SequenceNoiseResult:
    """Minimize the selected utility objective subject to finite exact MAP risk.

    For product channels, the full MAP success is the coordinate MAP success
    raised to coord_count, so the finite LP constrains the per-coordinate risk
    to target^(1 / coord_count).
    """
    noise = noise_support(kind, radius)
    index = {int(v): i for i, v in enumerate(noise)}
    ys = list(output_support(channel, noise))
    target_per_coord = math.exp(target_log2 * math.log(2.0) / coord_count)
    zero = np.zeros(len(noise), dtype=float)
    zero_idx = index.get(0)
    if zero_idx is not None:
        zero[zero_idx] = 1.0
        zero_risk = exact_map_risk(channel, noise, zero)
        if zero_risk <= target_per_coord:
            return _result(
                status="zero_noise_feasible",
                method="exact_map_lp",
                objective=0.0,
                support=noise,
                probabilities=zero,
                risk_per_coord=zero_risk,
                target_per_coord=target_per_coord,
                coord_count=coord_count,
                solver="none",
            )

    q = cp.Variable(len(noise), nonneg=True)
    t = cp.Variable(len(ys), nonneg=True)
    costs = utility_costs(noise, objective)
    risk_scale = float(np.max(channel.max_prior_by_leakage))
    if risk_scale <= 0.0:
        raise ValueError("channel max-prior values must have positive mass")
    constraints = [cp.sum(q) == 1.0]
    if kind == "symmetric":
        constraints.append(noise.astype(float) @ q == 0.0)

    for yi, y in enumerate(ys):
        for h, max_prior in zip(channel.leakage_values, channel.max_prior_by_leakage):
            j = index.get(int(y - h))
            if j is not None and max_prior > 0.0:
                constraints.append(t[yi] >= (float(max_prior) / risk_scale) * q[j])
    constraints.append(cp.sum(t) <= target_per_coord / risk_scale)

    problem = cp.Problem(cp.Minimize(costs @ q), constraints)
    problem.solve(solver=solver)
    probs = _clean_probabilities(q.value if q.value is not None else np.zeros(len(noise)))
    risk = exact_map_risk(channel, noise, probs)
    objective = float(np.sum(costs * probs))
    return _result(
        status=str(problem.status),
        method="exact_map_lp",
        objective=objective,
        support=noise,
        probabilities=probs,
        risk_per_coord=risk,
        target_per_coord=target_per_coord,
        coord_count=coord_count,
        solver=solver,
    )


def solve_exact_map_noise_scipy(
    channel: Channel,
    coord_count: int,
    target_log2: float,
    kind: str,
    radius: int,
    *,
    solver: str = "highs",
    objective: str = "second_moment",
) -> SequenceNoiseResult:
    """Sparse HiGHS implementation of the exact MAP LP.

    This solves the same LP as solve_exact_map_noise, but avoids CVXPY
    canonicalization overhead.  The MAP constraints are scaled by the largest
    prior atom to avoid numerical issues for 200+ bit secrets.
    """
    noise = noise_support(kind, radius)
    index = {int(v): i for i, v in enumerate(noise)}
    ys = list(output_support(channel, noise))
    n_q = len(noise)
    n_t = len(ys)
    n_var = n_q + n_t
    target_per_coord = math.exp(target_log2 * math.log(2.0) / coord_count)

    zero = np.zeros(n_q, dtype=float)
    zero_idx = index.get(0)
    if zero_idx is not None:
        zero[zero_idx] = 1.0
        zero_risk = exact_map_risk(channel, noise, zero)
        if zero_risk <= target_per_coord:
            return _result(
                status="zero_noise_feasible",
                method="exact_map_lp_scipy",
                objective=0.0,
                support=noise,
                probabilities=zero,
                risk_per_coord=zero_risk,
                target_per_coord=target_per_coord,
                coord_count=coord_count,
                solver="none",
            )

    risk_scale = float(np.max(channel.max_prior_by_leakage))
    if risk_scale <= 0.0:
        raise ValueError("channel max-prior values must have positive mass")

    rows = []
    rhs = []

    # sum_y t_y <= target / scale.
    row = {n_q + yi: 1.0 for yi in range(n_t)}
    rows.append(row)
    rhs.append(target_per_coord / risk_scale)

    # scaled MAP constraints: (m_h/scale) q_j - t_y <= 0.
    for yi, y in enumerate(ys):
        for h, max_prior in zip(channel.leakage_values, channel.max_prior_by_leakage):
            j = index.get(int(y - h))
            if j is not None and max_prior > 0.0:
                rows.append({j: float(max_prior) / risk_scale, n_q + yi: -1.0})
                rhs.append(0.0)

    a_ub = lil_matrix((len(rows), n_var), dtype=float)
    for i, row_entries in enumerate(rows):
        for j, value in row_entries.items():
            a_ub[i, j] = value

    a_eq_rows = 1 if kind != "symmetric" else 2
    a_eq = lil_matrix((a_eq_rows, n_var), dtype=float)
    a_eq[0, :n_q] = np.ones(n_q)
    b_eq = [1.0]
    if kind == "symmetric":
        a_eq[1, :n_q] = noise.astype(float)
        b_eq.append(0.0)

    costs = utility_costs(noise, objective)
    c = np.zeros(n_var, dtype=float)
    c[:n_q] = costs
    bounds = [(0.0, None)] * n_var

    opt = linprog(
        c,
        A_ub=a_ub.tocsr(),
        b_ub=np.asarray(rhs, dtype=float),
        A_eq=a_eq.tocsr(),
        b_eq=np.asarray(b_eq, dtype=float),
        bounds=bounds,
        method=solver,
    )
    probs = _clean_probabilities(opt.x[:n_q] if opt.x is not None else np.zeros(n_q))
    risk = exact_map_risk(channel, noise, probs)
    objective_value = float(np.sum(costs * probs))
    return _result(
        status=str(opt.message if not opt.success else "optimal"),
        method="exact_map_lp_scipy",
        objective=objective_value,
        support=noise,
        probabilities=probs,
        risk_per_coord=risk,
        target_per_coord=target_per_coord,
        coord_count=coord_count,
        solver=solver,
    )


def solve_prw_alpha_noise(
    channel: Channel,
    coord_count: int,
    target_log2: float,
    kind: str,
    radius: int,
    alpha: float,
    *,
    solver: str = "CLARABEL",
    objective: str = "second_moment",
) -> SequenceNoiseResult:
    """Minimize the selected utility objective subject to finite-alpha PRW."""
    if alpha < 1.0:
        raise ValueError("alpha must be at least one")

    noise = noise_support(kind, radius)
    index = {int(v): i for i, v in enumerate(noise)}
    ys = list(output_support(channel, noise))
    log_target_per_coord = target_log2 * math.log(2.0) / coord_count
    target_per_coord = math.exp(log_target_per_coord)
    zero = np.zeros(len(noise), dtype=float)
    zero_idx = index.get(0)
    if zero_idx is not None:
        zero[zero_idx] = 1.0
        zero_log_risk = log_prw_alpha_value(channel, noise, zero, alpha)
        if zero_log_risk <= log_target_per_coord:
            return _result(
                status="zero_noise_feasible",
                method=f"prw_alpha_{alpha:g}",
                objective=0.0,
                support=noise,
                probabilities=zero,
                risk_per_coord=math.exp(zero_log_risk) if zero_log_risk > -745.0 else 0.0,
                target_per_coord=target_per_coord,
                coord_count=coord_count,
                solver="none",
            )

    q = cp.Variable(len(noise), nonneg=True)
    costs = utility_costs(noise, objective)
    constraints = [cp.sum(q) == 1.0]
    if kind == "symmetric":
        constraints.append(noise.astype(float) @ q == 0.0)

    log_s_alpha = _log_sum_prior_alpha(channel, alpha)
    finite_log_s = log_s_alpha[np.isfinite(log_s_alpha)]
    if finite_log_s.size == 0:
        raise ValueError("channel alpha-prior weights must have positive mass")
    log_scale = float(np.max(finite_log_s))
    scaled_s_alpha = np.zeros_like(log_s_alpha, dtype=float)
    finite = np.isfinite(log_s_alpha)
    scaled_s_alpha[finite] = np.exp(log_s_alpha[finite] - log_scale)
    log_target_scaled = log_target_per_coord - log_scale / alpha
    if log_target_scaled > 700.0:
        target_scaled = math.exp(700.0)
    elif log_target_scaled < -745.0:
        target_scaled = 0.0
    else:
        target_scaled = math.exp(log_target_scaled)
    terms = []
    for y in ys:
        pieces = []
        for h, weight in zip(channel.leakage_values, scaled_s_alpha):
            j = index.get(int(y - h))
            if j is not None and weight > 0.0:
                pieces.append((float(weight) ** (1.0 / alpha)) * q[j])
        if pieces:
            terms.append(cp.norm(cp.hstack(pieces), p=alpha))
    constraints.append(cp.sum(cp.hstack(terms)) <= target_scaled)

    problem = cp.Problem(cp.Minimize(costs @ q), constraints)
    problem.solve(solver=solver)
    probs = _clean_probabilities(q.value if q.value is not None else np.zeros(len(noise)))
    risk = prw_alpha_value(channel, noise, probs, alpha)
    objective = float(np.sum(costs * probs))
    return _result(
        status=str(problem.status),
        method=f"prw_alpha_{alpha:g}",
        objective=objective,
        support=noise,
        probabilities=probs,
        risk_per_coord=risk,
        target_per_coord=target_per_coord,
        coord_count=coord_count,
        solver=solver,
    )


def summarize_distribution(result: SequenceNoiseResult, top: int = 12) -> str:
    order = np.argsort(-result.probabilities)
    parts = []
    for idx in order[:top]:
        p = float(result.probabilities[idx])
        if p <= 0.0:
            continue
        parts.append(f"{int(result.support[idx])}:{p:.4g}")
    return ", ".join(parts)


def _demo() -> None:
    rsa = make_rsa_timing_bit_channel(unit_ns=1.0)
    aes = make_aes_power_channel()
    cases = [
        ("RSA exact zero-mean", solve_exact_map_noise(rsa, 256, -150, "symmetric", 64)),
        ("RSA exact delay-only", solve_exact_map_noise(rsa, 256, -150, "onesided", 64)),
        ("AES exact zero-mean", solve_exact_map_noise(aes, 32, -150, "symmetric", 64)),
    ]
    for name, res in cases:
        print(name)
        print(
            f"  status={res.status} objective={res.objective:.6g} "
            f"mean={res.mean:.6g} risk={res.risk_per_coord:.6g} "
            f"target={res.target_per_coord:.6g}"
        )
        print(f"  top mass: {summarize_distribution(res)}")


if __name__ == "__main__":
    _demo()
