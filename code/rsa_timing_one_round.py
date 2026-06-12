"""One-round RSA timing-leakage sweep for empirical PRW calibration.

The secret is an l-bit random string.  We sample a finite support of such
strings, evaluate a deterministic RSA-style timing oracle on that support, and
then compute the Gaussian output-noise scale needed to certify several target
posterior success rates for different alpha values.

This script follows the first-round simplification requested for simulation:
sampled datapoints define the empirical expectation exactly; no validation set
or Hoeffding correction is used.  The PRW prior mass inside pi(x) is still the
true uniform l-bit prior, pi(x)=2^{-l}, so targets near 2^{-l} are meaningful.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from empirical_composition import EmpiricalAdaptiveComposer, as_2d_array


DEFAULT_ALPHAS = [16, 32, 64, 128, 256, 512, 1024]
DEFAULT_OFFSETS = [1, 2, 4, 8]


def rsa_timing_proxy(secret: Any, query: Any, t: int, transcript: Sequence[Any]) -> np.ndarray:
    """Deterministic RSA-style timing proxy for square-and-multiply.

    The query supplies a modulus and base.  The secret bit string is interpreted
    as the private exponent.  The returned scalar is a synthetic operation-cost
    trace: every bit incurs a square operation; one-bits incur an additional
    multiply operation; small operand-dependent terms make the leakage depend on
    the query/base as a black-box timing channel would.
    """

    bits = np.asarray(secret, dtype=np.uint8).reshape(-1)
    modulus = int(query["modulus"])
    base = int(query["base"]) % modulus
    mask = (1 << int(query["mask_bits"])) - 1

    acc = 1
    cost = 0.0
    for bit in bits:
        before = acc
        acc = (acc * acc) % modulus
        cost += 1.0 + 0.0005 * ((before ^ acc) & mask).bit_count()

        if bit:
            before = acc
            acc = (acc * base) % modulus
            cost += 1.2 + 0.0007 * ((before ^ acc ^ base) & mask).bit_count()
        else:
            cost += 0.05

    return np.array([cost], dtype=float)


def make_query(rng: np.random.Generator) -> dict[str, int]:
    """Create one fixed RSA-like query."""

    modulus = (1 << 257) - 93
    base = int(rng.integers(2, 1 << 63))
    return {"modulus": modulus, "base": base, "mask_bits": 64}


def target_log_r(secret_bits: int, alpha: float, offset: int) -> float:
    """Return log r for target delta=2^{-secret_bits+offset}.

    For one round under a true uniform l-bit prior:

        delta^alpha <= 2^{-l(alpha-1)} r.

    Therefore log_2 r = alpha * offset - l.
    """

    return (alpha * offset - secret_bits) * math.log(2.0)


def run_sweep(
    secret_bits: int,
    sample_count: int,
    alphas: Sequence[float],
    offsets: Sequence[int],
    seed: int,
    output_csv: Path,
) -> list[dict[str, str]]:
    rng = np.random.default_rng(seed)
    samples = rng.integers(0, 2, size=(sample_count, secret_bits), dtype=np.uint8)
    query = make_query(rng)
    means = as_2d_array(rsa_timing_proxy(sample, query, 1, []) for sample in samples)

    weights = np.full(sample_count, 1.0 / sample_count)
    prior_log_probs = np.full(sample_count, -secret_bits * math.log(2.0))

    leakage_mean = float(np.mean(means[:, 0]))
    leakage_std = float(np.std(means[:, 0]))
    leakage_min = float(np.min(means[:, 0]))
    leakage_max = float(np.max(means[:, 0]))

    rows: list[dict[str, str]] = []
    for alpha in alphas:
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

        no_leakage_log2_delta = -secret_bits * (alpha - 1.0) / alpha
        for offset in offsets:
            log_r = target_log_r(secret_bits, alpha, offset)
            target_log2_delta = -secret_bits + offset
            log2_r = log_r / math.log(2.0)

            row = {
                "secret_bits": str(secret_bits),
                "sample_count": str(sample_count),
                "alpha": f"{alpha:g}",
                "target": f"2^({target_log2_delta})",
                "target_log2_delta": f"{target_log2_delta:g}",
                "offset": str(offset),
                "no_leakage_log2_delta_bound": f"{no_leakage_log2_delta:.6f}",
                "log2_r": f"{log2_r:.6f}",
                "sigma": "",
                "sigma_over_leakage_std": "",
                "achieved_log2_delta_bound": "",
                "achieved_log2_local_ratio": "",
                "status": "",
            }

            if log_r < 0:
                row["status"] = "infeasible_target_below_no_leakage_alpha_bound"
            elif log_r == 0:
                row["status"] = "requires_infinite_noise_at_alpha_threshold"
            else:
                calibration = composer.calibrate_round(means, log_r_t=log_r)
                achieved_log2_delta = (
                    -secret_bits * (alpha - 1.0) + log2_r
                ) / alpha
                row["sigma"] = f"{calibration.sigma:.10f}"
                row["sigma_over_leakage_std"] = f"{calibration.sigma / leakage_std:.10f}"
                row["achieved_log2_delta_bound"] = f"{achieved_log2_delta:.6f}"
                row["achieved_log2_local_ratio"] = (
                    f"{calibration.log_local_ratio / math.log(2.0):.6f}"
                )
                row["status"] = "ok"

            rows.append(row)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    metadata = output_csv.with_suffix(".meta.txt")
    metadata.write_text(
        "\n".join(
            [
                f"secret_bits={secret_bits}",
                f"sample_count={sample_count}",
                f"seed={seed}",
                f"query_base={query['base']}",
                f"query_modulus={query['modulus']}",
                f"leakage_mean={leakage_mean:.10f}",
                f"leakage_std={leakage_std:.10f}",
                f"leakage_min={leakage_min:.10f}",
                f"leakage_max={leakage_max:.10f}",
            ]
        )
        + "\n"
    )
    return rows


def print_rows(rows: Sequence[dict[str, str]]) -> None:
    print(
        "alpha  target       log2_r      sigma        sigma/std    "
        "achieved_log2_delta  status"
    )
    for row in rows:
        sigma = row["sigma"] if row["sigma"] else "-"
        sigma_std = row["sigma_over_leakage_std"] if row["sigma_over_leakage_std"] else "-"
        achieved = row["achieved_log2_delta_bound"] if row["achieved_log2_delta_bound"] else "-"
        print(
            f"{row['alpha']:>5}  "
            f"{row['target']:>11}  "
            f"{row['log2_r']:>10}  "
            f"{sigma:>11}  "
            f"{sigma_std:>10}  "
            f"{achieved:>20}  "
            f"{row['status']}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secret-bits", type=int, default=256)
    parser.add_argument("--sample-count", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20260603)
    parser.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        default=DEFAULT_ALPHAS,
        help="alpha values to sweep",
    )
    parser.add_argument(
        "--offsets",
        type=int,
        nargs="+",
        default=DEFAULT_OFFSETS,
        help="targets are 2^(-secret_bits + offset)",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "rsa_timing_one_round.csv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = run_sweep(
        secret_bits=args.secret_bits,
        sample_count=args.sample_count,
        alphas=args.alphas,
        offsets=args.offsets,
        seed=args.seed,
        output_csv=args.output_csv,
    )
    print_rows(rows)
    print(f"\nwrote {args.output_csv}")
    print(f"wrote {args.output_csv.with_suffix('.meta.txt')}")


if __name__ == "__main__":
    main()
