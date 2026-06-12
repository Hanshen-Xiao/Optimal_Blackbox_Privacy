"""Estimate validation data needs for confidence-aware RSA timing simulation.

This script answers: after the no-confidence empirical calibration chooses a
Gaussian scale, how many independent validation samples would Hoeffding need to
certify the local condition with confidence parameter gamma?

The calculation uses a target validation margin equal to a user-selected
fraction of the local budget:

    beta <= margin_fraction * r_t.

For large alpha the required sample count can be astronomically large, so the
script reports log10(m_validation) instead of only integer counts.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from confidence_composition import log10_required_validation_samples
from empirical_composition import EmpiricalAdaptiveComposer, as_2d_array, gaussian_logpdf
from rsa_timing_one_round import make_query, rsa_timing_proxy


def log_sub_exp(log_a: float, log_b: float) -> float:
    """Return log(exp(log_a)-exp(log_b)) for log_a >= log_b."""

    if log_b > log_a:
        raise ValueError("log_b must be <= log_a")
    if log_a == log_b:
        return float("-inf")
    return log_a + math.log1p(-math.exp(log_b - log_a))


def log_range_from_log_values(log_values: np.ndarray) -> float:
    """Log range for positive values represented in log-space."""

    max_log = float(np.max(log_values))
    min_log = float(np.min(log_values))
    if max_log == min_log:
        return float("-inf")
    return log_sub_exp(max_log, min_log)


def precompute(samples: np.ndarray, queries: Sequence[dict[str, int]]) -> list[np.ndarray]:
    return [
        as_2d_array(rsa_timing_proxy(sample, query, t + 1, []) for sample in samples)
        for t, query in enumerate(queries)
    ]


def run_requirement_sweep(
    secret_bits: int,
    sample_count: int,
    alphas: Sequence[float],
    offset: int,
    rounds: int,
    gamma: float,
    max_candidates: int,
    margin_fraction: float,
    seed: int,
    output_csv: Path,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    samples = rng.integers(0, 2, size=(sample_count, secret_bits), dtype=np.uint8)
    queries = [make_query(rng) for _ in range(rounds)]
    means_by_round = precompute(samples, queries)
    actual_secret = samples[int(rng.integers(sample_count))]

    rows: list[dict[str, object]] = []
    weights = np.full(sample_count, 1.0 / sample_count)
    prior_log_probs = np.full(sample_count, -secret_bits * math.log(2.0))
    all_means = np.vstack(means_by_round)[:, 0]
    leakage_std = float(np.std(all_means))

    for alpha in alphas:
        total_log2_r = alpha * offset - secret_bits
        row_prefix = {
            "secret_bits": secret_bits,
            "sample_count_no_confidence": sample_count,
            "alpha": f"{alpha:g}",
            "offset": offset,
            "rounds": rounds,
            "gamma": gamma,
            "max_candidates": max_candidates,
            "margin_fraction": margin_fraction,
            "log2_product_r": total_log2_r,
            "log2_r_per_round": total_log2_r / rounds,
        }
        if total_log2_r <= 0:
            rows.append(
                {
                    **row_prefix,
                    "status": "infeasible_or_infinite_noise_threshold",
                    "max_log10_validation_samples_per_round": "",
                    "total_log10_blackbox_evals": "",
                    "largest_round": "",
                }
            )
            continue

        log_r_per_round = (total_log2_r / rounds) * math.log(2.0)
        composer = EmpiricalAdaptiveComposer(
            samples=samples,
            leakage_fn=rsa_timing_proxy,
            alpha=alpha,
            weights=weights,
            prior_log_probs=prior_log_probs,
            rng=rng,
            sigma_initial_hi=max(1.0, leakage_std),
            sigma_tol=1e-7,
        )

        per_round_log10_m: list[float] = []
        sigmas: list[float] = []
        for t, (query, means) in enumerate(zip(queries, means_by_round), start=1):
            calibration = composer.calibrate_round(means, log_r_t=log_r_per_round)
            sigmas.append(calibration.sigma)

            sq_dist = np.sum((means - calibration.reference_mean.reshape(1, -1)) ** 2, axis=1)
            log_div = alpha * (alpha - 1.0) * sq_dist / (2.0 * calibration.sigma**2)
            log_range = log_range_from_log_values(log_div)
            log_margin = math.log(margin_fraction) + log_r_per_round
            log10_m = log10_required_validation_samples(
                log_range,
                log_margin,
                gamma,
                max_candidates,
            )
            per_round_log10_m.append(log10_m)

            raw_output = rsa_timing_proxy(actual_secret, query, t, [])
            noisy_output = raw_output + rng.normal(0.0, calibration.sigma, size=raw_output.shape)
            composer._update_prefix(  # pylint: disable=protected-access
                means,
                calibration.reference_mean,
                calibration.sigma,
                noisy_output,
            )
            composer.accountant.add_log_round(log_r_per_round)

        max_log10_m = max(per_round_log10_m)
        largest_round = int(np.argmax(per_round_log10_m) + 1)
        # If every round uses this validation size plus one equally-sized
        # training set, total leakage evaluations scale as 2*T*m.
        total_log10 = math.log10(2.0 * rounds) + max_log10_m
        rows.append(
            {
                **row_prefix,
                "status": "ok",
                "max_log10_validation_samples_per_round": f"{max_log10_m:.3f}",
                "total_log10_blackbox_evals": f"{total_log10:.3f}",
                "largest_round": largest_round,
                "avg_sigma2": f"{float(np.mean(np.square(sigmas))):.6f}",
            }
        )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def print_rows(rows: Sequence[dict[str, object]]) -> None:
    print("alpha  log2_r/t  log10(m_val/round)  log10(total evals)  status")
    for row in rows:
        print(
            f"{row['alpha']:>5}  "
            f"{float(row['log2_r_per_round']):>8.1f}  "
            f"{str(row['max_log10_validation_samples_per_round']):>20}  "
            f"{str(row['total_log10_blackbox_evals']):>18}  "
            f"{row['status']}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secret-bits", type=int, default=256)
    parser.add_argument("--sample-count", type=int, default=4096)
    parser.add_argument("--alphas", type=float, nargs="+", default=[128, 256, 512, 1024, 2048, 4096])
    parser.add_argument("--offset", type=int, default=8)
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--gamma", type=float, default=0.05)
    parser.add_argument("--max-candidates", type=int, default=20)
    parser.add_argument("--margin-fraction", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=20260603)
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "rsa_confidence_sample_requirements.csv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = run_requirement_sweep(
        secret_bits=args.secret_bits,
        sample_count=args.sample_count,
        alphas=args.alphas,
        offset=args.offset,
        rounds=args.rounds,
        gamma=args.gamma,
        max_candidates=args.max_candidates,
        margin_fraction=args.margin_fraction,
        seed=args.seed,
        output_csv=args.output_csv,
    )
    print_rows(rows)
    print(f"\nwrote {args.output_csv}")


if __name__ == "__main__":
    main()
