"""Short-secret aggregate-power experiment with PAC comparison."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from convex_replot_experiments import line_plot, write_csv
from discrete_leakage_experiments import binary_kl, mutual_information_coord
from discrete_noise_convex import (
    make_binomial_aggregate_channel,
    noise_support,
    solve_exact_map_noise_scipy,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "power_short"


def discrete_gaussian_on_support(support: np.ndarray, scale: float) -> np.ndarray:
    if scale <= 0.0:
        q = np.zeros(len(support), dtype=float)
        zero_idx = np.where(support == 0)[0]
        if len(zero_idx) == 0:
            raise ValueError("support must contain zero")
        q[int(zero_idx[0])] = 1.0
        return q
    logq = -(support.astype(float) ** 2) / (2.0 * scale * scale)
    logq -= float(np.max(logq))
    q = np.exp(logq)
    q /= float(np.sum(q))
    return q


def second_moment(support: np.ndarray, q: np.ndarray) -> float:
    return float(np.sum((support.astype(float) ** 2) * q))


def mean_noise(support: np.ndarray, q: np.ndarray) -> float:
    return float(np.sum(support.astype(float) * q))


def calibrate_pac_gaussian_unconstrained(channel, target_log2: float) -> dict:
    """PAC/MI calibration for unconstrained continuous Gaussian noise.

    For H distributed as the aggregate leakage and Y = H + N with
    N ~ Gaussian(0, sigma^2), use the AWGN capacity bound

        I(H;Y) <= 0.5 log(1 + Var(H) / sigma^2).

    This is the standard PAC-style MI comparison without finite-support,
    zero-mean-sequence, or one-sided constraints.
    """
    prior = channel.prior_success_coord
    target = 2.0 ** target_log2
    threshold = binary_kl(target, prior)
    values = channel.leakage_values.astype(float)
    mean = float(np.sum(channel.leakage_probs * values))
    variance = float(np.sum(channel.leakage_probs * (values - mean) ** 2))
    if threshold <= 0.0:
        return {
            "status": "infeasible_threshold",
            "second_moment": float("inf"),
            "rms": float("inf"),
            "mean": 0.0,
            "scale": float("inf"),
            "mi_bound_nats": "",
            "mi_threshold_nats": threshold,
        }
    sigma2 = variance / math.expm1(2.0 * threshold)
    mi_bound = 0.5 * math.log1p(variance / sigma2)
    return {
        "status": "ok",
        "second_moment": sigma2,
        "rms": math.sqrt(sigma2),
        "mean": 0.0,
        "scale": math.sqrt(sigma2),
        "mi_bound_nats": mi_bound,
        "mi_threshold_nats": threshold,
    }


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lengths = [128, 127, 126, 125, 124]
    targets = [-123.0, -122.0]
    supports = {
        "zero_mean": ("symmetric", 75, noise_support("symmetric", 75)),
        "less_power_only": ("nonpositive", 150, noise_support("nonpositive", 150)),
    }
    rows: List[dict] = []

    for bits in lengths:
        channel = make_binomial_aggregate_channel(bits, name=f"AES aggregate HW({bits})")
        no_noise_log2 = math.log2((bits + 1) * 2.0 ** (-bits))
        for target in targets:
            for constraint, (kind, radius, support) in supports.items():
                exact = solve_exact_map_noise_scipy(channel, 1, target, kind, radius)
                rows.append({
                    "experiment": "power_short",
                    "method": "Exact MAP sequence LP",
                    "secret_bits": bits,
                    "target_log2": target,
                    "constraint": constraint,
                    "status": exact.status,
                    "second_moment": exact.objective,
                    "rms_power_bins": math.sqrt(max(exact.objective, 0.0)),
                    "mean_noise": exact.mean,
                    "risk_log2": exact.full_log_risk / math.log(2.0),
                    "no_noise_log2": no_noise_log2,
                    "prior_log2": -float(bits),
                    "support_size": len(support),
                    "support_min": int(support[0]),
                    "support_max": int(support[-1]),
                    "mi_nats": "",
                    "mi_threshold_nats": "",
                })

            pac = calibrate_pac_gaussian_unconstrained(channel, target)
            rows.append({
                "experiment": "power_short",
                "method": "PAC MI unconstrained Gaussian",
                "secret_bits": bits,
                "target_log2": target,
                "constraint": "unconstrained",
                "status": pac["status"],
                "second_moment": pac["second_moment"],
                "rms_power_bins": pac["rms"],
                "mean_noise": pac["mean"],
                "risk_log2": "",
                "no_noise_log2": no_noise_log2,
                "prior_log2": -float(bits),
                "support_size": "continuous",
                "support_min": "-inf",
                "support_max": "inf",
                "mi_nats": pac["mi_bound_nats"],
                "mi_threshold_nats": pac["mi_threshold_nats"],
            })

    write_csv(OUT / "power_short_results.csv", rows)

    exact_series: Dict[str, List[Tuple[float, float]]] = {}
    for target in targets:
        for constraint, label in [
            ("zero_mean", "zero-mean"),
            ("less_power_only", "less-only"),
        ]:
            key = f"{label}, 2^{int(target)}"
            exact_series[key] = sorted([
                (float(r["secret_bits"]), float(r["second_moment"]))
                for r in rows
                if r["method"] == "Exact MAP sequence LP"
                and r["constraint"] == constraint
                and r["target_log2"] == target
            ])

    line_plot(
        OUT / "power_short_exact_noise.png",
        "Aggregate AES Power: Exact Noise Optimization",
        "151-variable finite sequence; scalar Hamming-weight leakage",
        "Secret length l (bits)",
        "Noise second moment (power-bin^2)",
        exact_series,
        log_y=True,
        x_ticks_override=lengths,
    )

    for target in targets:
        compare_series: Dict[str, List[Tuple[float, float]]] = {}
        for constraint, label in [
            ("zero_mean", "zero-mean"),
            ("less_power_only", "less-only"),
        ]:
            compare_series[f"{label} exact"] = sorted([
                (float(r["secret_bits"]), math.log10(float(r["second_moment"])))
                for r in rows
                if r["method"] == "Exact MAP sequence LP"
                and r["constraint"] == constraint
                and r["target_log2"] == target
            ])
        compare_series["PAC Gaussian"] = sorted([
            (float(r["secret_bits"]), math.log10(float(r["second_moment"])))
            for r in rows
            if r["method"] == "PAC MI unconstrained Gaussian"
            and r["target_log2"] == target
        ])

        line_plot(
            OUT / f"power_short_pac_comparison_target{abs(int(target))}.png",
            f"Aggregate AES Power: PAC Comparison for Target 2^{int(target)}",
            "PAC uses unconstrained continuous Gaussian noise and an MI bound",
            "Secret length l (bits)",
            "log10 noise second moment",
            compare_series,
            log_y=False,
            x_ticks_override=lengths,
        )

    exponent_series: Dict[str, List[Tuple[float, float]]] = {
        "no-noise exponent": [],
        "target 2^-123": [],
        "target 2^-122": [],
        "prior exponent": [],
    }
    seen = set()
    for row in rows:
        bits = int(row["secret_bits"])
        if bits in seen:
            continue
        seen.add(bits)
        exponent_series["no-noise exponent"].append((bits, -float(row["no_noise_log2"])))
        exponent_series["prior exponent"].append((bits, -float(row["prior_log2"])))
        exponent_series["target 2^-123"].append((bits, 123.0))
        exponent_series["target 2^-122"].append((bits, 122.0))

    line_plot(
        OUT / "power_short_posterior_exponents.png",
        "Aggregate AES Power: Posterior Success Exponents",
        "Higher exponent means lower posterior success",
        "Secret length l (bits)",
        "-log2 posterior success",
        {k: sorted(v) for k, v in exponent_series.items()},
        x_ticks_override=lengths,
        y_min_override=116.0,
        y_max_override=130.0,
    )

    print(f"wrote short power results to {OUT}")


if __name__ == "__main__":
    run()
