"""Non-uniform aggregate-power experiment with closed-form ground truth."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from convex_replot_experiments import line_plot, write_csv
from discrete_leakage_experiments import binary_kl
from discrete_noise_convex import (
    make_binomial_aggregate_channel,
    noise_support,
    solve_exact_map_noise_scipy,
)
from power_short_experiments import calibrate_pac_gaussian_unconstrained


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "power_nonuniform"


def logD_alpha_bernoulli(rho: float, alpha: float, q: float) -> float:
    eps = 1e-300
    rho = min(max(float(rho), eps), 1.0 - eps)
    q = min(max(float(q), eps), 1.0 - eps)
    log_a = math.log(q) + alpha * (math.log(rho) - math.log(q))
    log_b = math.log1p(-q) + alpha * (math.log1p(-rho) - math.log1p(-q))
    return float(logsumexp([log_a, log_b]))


def pac_alpha_log_cost(bits: int, p: float, sigma2: float, alpha: float) -> float:
    channel = make_binomial_aggregate_channel(bits, bit_prob_one=p)
    h = channel.leakage_values.astype(float)
    mu = bits * p
    c = alpha * (alpha - 1.0) / 2.0
    log_probs = np.log(channel.leakage_probs)
    return float(logsumexp(log_probs + c * ((h - mu) ** 2) / sigma2))


def pac_alpha_sigma2(bits: int, p: float, target_log2: float, alpha: float) -> Tuple[float, float]:
    channel = make_binomial_aggregate_channel(bits, bit_prob_one=p)
    rho = 2.0 ** target_log2
    q = channel.prior_success_coord
    target = logD_alpha_bernoulli(rho, alpha, q)
    if target <= 0.0:
        return float("inf"), target
    lo, hi = 1e-12, 1.0
    while pac_alpha_log_cost(bits, p, hi, alpha) > target:
        hi *= 2.0
        if hi > 1e30:
            return float("inf"), target
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pac_alpha_log_cost(bits, p, mid, alpha) <= target:
            hi = mid
        else:
            lo = mid
    return hi, target


def optimize_pac_alpha(bits: int, p: float, target_log2: float) -> dict:
    channel = make_binomial_aggregate_channel(bits, bit_prob_one=p)
    rho = 2.0 ** target_log2
    if rho <= channel.prior_success_coord:
        return {
            "status": "target_below_prior",
            "alpha": "",
            "second_moment": float("inf"),
            "rms": float("inf"),
            "threshold": "",
        }

    def objective(log_alpha_minus_one: float) -> float:
        alpha = 1.0 + math.exp(log_alpha_minus_one)
        sigma2, _ = pac_alpha_sigma2(bits, p, target_log2, alpha)
        if not math.isfinite(sigma2) or sigma2 <= 0.0:
            return 1e300
        return math.log(sigma2)

    best = None
    # This broad range covers alpha close to one through very large moments.
    for lo, hi in [(-4.0, 4.0), (2.0, 7.0), (5.0, 10.0)]:
        res = minimize_scalar(objective, bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
        if best is None or res.fun < best.fun:
            best = res
    alpha = 1.0 + math.exp(float(best.x))
    sigma2, threshold = pac_alpha_sigma2(bits, p, target_log2, alpha)
    return {
        "status": "ok",
        "alpha": alpha,
        "second_moment": sigma2,
        "rms": math.sqrt(sigma2),
        "threshold": threshold,
    }


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lengths = [128, 127, 126, 125, 124]
    targets = [-120.0, -119.0]
    p = 0.51
    supports = {
        "zero_mean": ("symmetric", 75, noise_support("symmetric", 75)),
        "less_power_only": ("nonpositive", 150, noise_support("nonpositive", 150)),
    }
    rows: List[dict] = []

    for bits in lengths:
        channel = make_binomial_aggregate_channel(bits, bit_prob_one=p)
        no_noise_log2 = math.log2(float(np.sum(channel.max_prior_by_leakage)))
        prior_log2 = math.log2(channel.prior_success_coord)
        shannon_entropy = -bits * (p * math.log2(p) + (1.0 - p) * math.log2(1.0 - p))
        for target in targets:
            for constraint, (kind, radius, support) in supports.items():
                exact = solve_exact_map_noise_scipy(channel, 1, target, kind, radius)
                rows.append({
                    "experiment": "power_nonuniform",
                    "method": "Exact MAP sequence LP",
                    "secret_bits": bits,
                    "bit_prob_one": p,
                    "shannon_entropy_bits": shannon_entropy,
                    "target_log2": target,
                    "constraint": constraint,
                    "status": exact.status,
                    "second_moment": exact.objective,
                    "rms_power_bins": math.sqrt(max(exact.objective, 0.0)),
                    "mean_noise": exact.mean,
                    "risk_log2": exact.full_log_risk / math.log(2.0),
                    "no_noise_log2": no_noise_log2,
                    "prior_log2": prior_log2,
                    "support_size": len(support),
                    "sample_count": "closed_form",
                    "pac_alpha": "",
                    "pac_threshold": "",
                })

            pac_mi = calibrate_pac_gaussian_unconstrained(channel, target)
            rows.append({
                "experiment": "power_nonuniform",
                "method": "PAC MI unconstrained Gaussian",
                "secret_bits": bits,
                "bit_prob_one": p,
                "shannon_entropy_bits": shannon_entropy,
                "target_log2": target,
                "constraint": "unconstrained",
                "status": pac_mi["status"],
                "second_moment": pac_mi["second_moment"],
                "rms_power_bins": pac_mi["rms"],
                "mean_noise": 0.0,
                "risk_log2": "",
                "no_noise_log2": no_noise_log2,
                "prior_log2": prior_log2,
                "support_size": "continuous",
                "sample_count": "closed_form",
                "pac_alpha": "",
                "pac_threshold": pac_mi["mi_threshold_nats"],
            })

            pac_alpha = optimize_pac_alpha(bits, p, target)
            rows.append({
                "experiment": "power_nonuniform",
                "method": "PAC alpha unconstrained Gaussian",
                "secret_bits": bits,
                "bit_prob_one": p,
                "shannon_entropy_bits": shannon_entropy,
                "target_log2": target,
                "constraint": "unconstrained",
                "status": pac_alpha["status"],
                "second_moment": pac_alpha["second_moment"],
                "rms_power_bins": pac_alpha["rms"],
                "mean_noise": 0.0,
                "risk_log2": "",
                "no_noise_log2": no_noise_log2,
                "prior_log2": prior_log2,
                "support_size": "continuous",
                "sample_count": "closed_form",
                "pac_alpha": pac_alpha["alpha"],
                "pac_threshold": pac_alpha["threshold"],
            })

    write_csv(OUT / "power_nonuniform_results.csv", rows)

    exact_series: Dict[str, List[Tuple[float, float]]] = {}
    for target in targets:
        for constraint, label in [("zero_mean", "zero-mean"), ("less_power_only", "less-only")]:
            exact_series[f"{label}, 2^{int(target)}"] = sorted([
                (float(r["secret_bits"]), max(float(r["second_moment"]), 1e-9))
                for r in rows
                if r["method"] == "Exact MAP sequence LP"
                and r["constraint"] == constraint
                and r["target_log2"] == target
            ])
    line_plot(
        OUT / "power_nonuniform_exact_noise.png",
        "Non-Uniform Aggregate Power: Exact Noise Optimization",
        "Bernoulli(0.51) bits, closed-form ground truth",
        "Secret length l (bits)",
        "Noise second moment (power-bin^2)",
        exact_series,
        log_y=True,
        x_ticks_override=lengths,
    )

    for target in targets:
        compare: Dict[str, List[Tuple[float, float]]] = {}
        for method, label in [
            ("Exact MAP sequence LP", "exact zero-mean"),
            ("PAC alpha unconstrained Gaussian", "PAC alpha"),
            ("PAC MI unconstrained Gaussian", "PAC MI"),
        ]:
            compare[label] = sorted([
                (float(r["secret_bits"]), math.log10(max(float(r["second_moment"]), 1e-9)))
                for r in rows
                if r["method"] == method
                and r["target_log2"] == target
                and (r["constraint"] in {"zero_mean", "unconstrained"})
            ])
        line_plot(
            OUT / f"power_nonuniform_pac_target{abs(int(target))}.png",
            f"Non-Uniform Power: PAC Comparison for Target 2^{int(target)}",
            "PAC curves use unconstrained continuous Gaussian noise",
            "Secret length l (bits)",
            "log10 noise second moment",
            compare,
            x_ticks_override=lengths,
        )

    exponents: Dict[str, List[Tuple[float, float]]] = {
        "no-noise exponent": [],
        "target 2^-120": [],
        "target 2^-119": [],
        "prior exponent": [],
    }
    seen = set()
    for row in rows:
        bits = int(row["secret_bits"])
        if bits in seen:
            continue
        seen.add(bits)
        exponents["no-noise exponent"].append((bits, -float(row["no_noise_log2"])))
        exponents["prior exponent"].append((bits, -float(row["prior_log2"])))
        exponents["target 2^-120"].append((bits, 120.0))
        exponents["target 2^-119"].append((bits, 119.0))
    line_plot(
        OUT / "power_nonuniform_exponents.png",
        "Non-Uniform Aggregate Power: Closed-Form Ground Truth",
        "Prior and no-noise posterior success are computed exactly",
        "Secret length l (bits)",
        "-log2 success probability",
        {k: sorted(v) for k, v in exponents.items()},
        x_ticks_override=lengths,
        y_min_override=112.0,
        y_max_override=126.0,
    )

    print(f"wrote non-uniform power results to {OUT}")


if __name__ == "__main__":
    run()
