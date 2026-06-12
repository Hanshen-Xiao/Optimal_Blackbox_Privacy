"""RSA timing-leakage composition sweep.

For each target posterior success rate, alpha, and total composition length T,
this script splits the total PRW alpha budget evenly across rounds and reports
the final per-round average Gaussian noise variance:

    (1/T) sum_{t=1}^T sigma_t^2.

It uses the same first-simulation simplification as rsa_timing_one_round.py:
sampled secrets define the empirical expectation exactly, with no validation
set and no Hoeffding correction.  The PRW prior mass is the true uniform l-bit
prior, pi(x)=2^{-l}.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from empirical_composition import EmpiricalAdaptiveComposer, as_2d_array
from rsa_timing_one_round import DEFAULT_OFFSETS, make_query, rsa_timing_proxy


DEFAULT_ALPHAS = [128, 256, 512, 1024]
DEFAULT_T_VALUES = [1, 2, 4, 8]


def log2_product_r(secret_bits: int, alpha: float, offset: int) -> float:
    """For target delta=2^{-secret_bits+offset}, return log2 prod_t r_t."""

    return alpha * offset - secret_bits


def precompute_leakage(
    samples: np.ndarray,
    queries: Sequence[dict[str, int]],
) -> list[np.ndarray]:
    return [
        as_2d_array(rsa_timing_proxy(sample, query, t + 1, []) for sample in samples)
        for t, query in enumerate(queries)
    ]


def run_precomputed_composition(
    samples: np.ndarray,
    means_by_round: Sequence[np.ndarray],
    queries: Sequence[dict[str, int]],
    alpha: float,
    log_r_per_round: float,
    secret_bits: int,
    rng: np.random.Generator,
) -> list[dict[str, float]]:
    sample_count = len(samples)
    weights = np.full(sample_count, 1.0 / sample_count)
    prior_log_probs = np.full(sample_count, -secret_bits * math.log(2.0))
    leakage_std = float(np.std(np.vstack(means_by_round)[:, 0]))

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

    actual_index = int(rng.integers(sample_count))
    actual_secret = samples[actual_index]
    records: list[dict[str, float]] = []
    sigmas: list[float] = []

    for t, (query, means) in enumerate(zip(queries, means_by_round), start=1):
        calibration = composer.calibrate_round(means, log_r_t=log_r_per_round)
        raw_output = rsa_timing_proxy(actual_secret, query, t, [])
        noise = rng.normal(0.0, calibration.sigma, size=raw_output.shape)
        noisy_output = raw_output + noise

        composer._update_prefix(  # pylint: disable=protected-access
            means,
            calibration.reference_mean,
            calibration.sigma,
            noisy_output,
        )
        composer.accountant.add_log_round(log_r_per_round)
        sigmas.append(calibration.sigma)

        records.append(
            {
                "round": float(t),
                "sigma": float(calibration.sigma),
                "avg_sigma2": float(np.mean(np.square(sigmas))),
                "rms_sigma": float(math.sqrt(np.mean(np.square(sigmas)))),
                "log2_delta_bound": float(composer.accountant.log_delta_bound / math.log(2.0)),
                "log2_local_ratio": float(calibration.log_local_ratio / math.log(2.0)),
            }
        )

    return records


def summarize_trials(
    trial_records: Sequence[Sequence[dict[str, float]]],
) -> list[dict[str, float]]:
    round_count = len(trial_records[0])
    rows: list[dict[str, float]] = []
    for idx in range(round_count):
        keys = trial_records[0][idx].keys()
        row: dict[str, float] = {}
        for key in keys:
            values = np.array([records[idx][key] for records in trial_records], dtype=float)
            row[key] = float(np.mean(values))
            if key in {"sigma", "avg_sigma2", "rms_sigma"}:
                row[f"{key}_std"] = float(np.std(values))
        rows.append(row)
    return rows


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_sweep(
    secret_bits: int,
    sample_count: int,
    alphas: Sequence[float],
    offsets: Sequence[int],
    t_values: Sequence[int],
    trials: int,
    seed: int,
    output_dir: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, float]]:
    rng = np.random.default_rng(seed)
    max_rounds = max(t_values)
    samples = rng.integers(0, 2, size=(sample_count, secret_bits), dtype=np.uint8)
    queries = [make_query(rng) for _ in range(max_rounds)]
    means_by_round = precompute_leakage(samples, queries)

    all_means = np.vstack(means_by_round)[:, 0]
    meta = {
        "secret_bits": float(secret_bits),
        "sample_count": float(sample_count),
        "trials": float(trials),
        "seed": float(seed),
        "leakage_mean": float(np.mean(all_means)),
        "leakage_std": float(np.std(all_means)),
        "leakage_min": float(np.min(all_means)),
        "leakage_max": float(np.max(all_means)),
    }

    final_rows: list[dict[str, object]] = []
    trajectory_rows: list[dict[str, object]] = []

    for alpha in alphas:
        for offset in offsets:
            total_log2_r = log2_product_r(secret_bits, alpha, offset)
            target_log2_delta = -secret_bits + offset

            for total_rounds in t_values:
                row: dict[str, object] = {
                    "secret_bits": secret_bits,
                    "sample_count": sample_count,
                    "trials": trials,
                    "alpha": f"{alpha:g}",
                    "target": f"2^({target_log2_delta})",
                    "target_log2_delta": target_log2_delta,
                    "offset": offset,
                    "total_rounds": total_rounds,
                    "log2_product_r": total_log2_r,
                    "log2_r_per_round": total_log2_r / total_rounds,
                    "final_avg_sigma2": "",
                    "final_rms_sigma": "",
                    "final_avg_sigma2_std": "",
                    "final_log2_delta_bound": "",
                    "status": "",
                }

                if total_log2_r < 0:
                    row["status"] = "infeasible_target_below_no_leakage_alpha_bound"
                    final_rows.append(row)
                    continue
                if total_log2_r == 0:
                    row["status"] = "requires_infinite_noise_at_alpha_threshold"
                    final_rows.append(row)
                    continue

                trial_records = [
                    run_precomputed_composition(
                        samples=samples,
                        means_by_round=means_by_round[:total_rounds],
                        queries=queries[:total_rounds],
                        alpha=alpha,
                        log_r_per_round=(total_log2_r / total_rounds) * math.log(2.0),
                        secret_bits=secret_bits,
                        rng=rng,
                    )
                    for _ in range(trials)
                ]
                summary = summarize_trials(trial_records)
                final = summary[-1]
                row["final_avg_sigma2"] = f"{final['avg_sigma2']:.10f}"
                row["final_rms_sigma"] = f"{final['rms_sigma']:.10f}"
                row["final_avg_sigma2_std"] = f"{final.get('avg_sigma2_std', 0.0):.10f}"
                row["final_log2_delta_bound"] = f"{final['log2_delta_bound']:.6f}"
                row["status"] = "ok"
                final_rows.append(row)

                if total_rounds == max_rounds:
                    for round_row in summary:
                        trajectory_rows.append(
                            {
                                "secret_bits": secret_bits,
                                "sample_count": sample_count,
                                "trials": trials,
                                "alpha": f"{alpha:g}",
                                "target": f"2^({target_log2_delta})",
                                "offset": offset,
                                "total_rounds": total_rounds,
                                "round": int(round_row["round"]),
                                "sigma": f"{round_row['sigma']:.10f}",
                                "avg_sigma2": f"{round_row['avg_sigma2']:.10f}",
                                "rms_sigma": f"{round_row['rms_sigma']:.10f}",
                                "log2_delta_bound": f"{round_row['log2_delta_bound']:.6f}",
                                "log2_local_ratio": f"{round_row['log2_local_ratio']:.6f}",
                            }
                        )

    write_csv(output_dir / "rsa_timing_composition_by_T.csv", final_rows)
    if trajectory_rows:
        write_csv(output_dir / f"rsa_timing_composition_trajectory_T{max_rounds}.csv", trajectory_rows)
    (output_dir / "rsa_timing_composition.meta.txt").write_text(
        "\n".join(f"{key}={value}" for key, value in meta.items()) + "\n"
    )

    return final_rows, trajectory_rows, meta


def print_compact(rows: Sequence[dict[str, object]], offset: int) -> None:
    ok_rows = [row for row in rows if row["offset"] == offset]
    print(f"Compact view for target offset={offset}: delta=2^(-256+{offset})")
    print("alpha  T   log2_r/t   final_avg_sigma2   final_rms_sigma   status")
    for row in ok_rows:
        print(
            f"{row['alpha']:>5}  "
            f"{row['total_rounds']:>1}  "
            f"{float(row['log2_r_per_round']):>10.3f}  "
            f"{str(row['final_avg_sigma2']):>17}  "
            f"{str(row['final_rms_sigma']):>15}  "
            f"{row['status']}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secret-bits", type=int, default=256)
    parser.add_argument("--sample-count", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20260603)
    parser.add_argument("--trials", type=int, default=1)
    parser.add_argument("--alphas", type=float, nargs="+", default=DEFAULT_ALPHAS)
    parser.add_argument("--offsets", type=int, nargs="+", default=DEFAULT_OFFSETS)
    parser.add_argument("--t-values", type=int, nargs="+", default=DEFAULT_T_VALUES)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs",
    )
    parser.add_argument("--compact-offset", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    final_rows, _, _ = run_sweep(
        secret_bits=args.secret_bits,
        sample_count=args.sample_count,
        alphas=args.alphas,
        offsets=args.offsets,
        t_values=args.t_values,
        trials=args.trials,
        seed=args.seed,
        output_dir=args.output_dir,
    )
    print_compact(final_rows, args.compact_offset)
    print(f"\nwrote {args.output_dir / 'rsa_timing_composition_by_T.csv'}")
    print(f"wrote {args.output_dir / ('rsa_timing_composition_trajectory_T' + str(max(args.t_values)) + '.csv')}")
    print(f"wrote {args.output_dir / 'rsa_timing_composition.meta.txt'}")


if __name__ == "__main__":
    main()
