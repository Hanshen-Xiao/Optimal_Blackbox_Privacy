"""Focused aggregate-power experiment near the scalar leakage boundary."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Tuple

from convex_replot_experiments import line_plot, write_csv
from discrete_noise_convex import (
    make_binomial_aggregate_channel,
    solve_exact_map_noise_scipy,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "power_focus"


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lengths = [248, 252, 254, 256]
    target_log2 = -247.0
    radius = 256
    rows: List[dict] = []

    for bits in lengths:
        channel = make_binomial_aggregate_channel(bits, name=f"AES aggregate HW({bits})")
        no_noise_log2 = math.log2((bits + 1) * 2.0 ** (-bits))
        prior_log2 = -float(bits)

        for kind, label in [
            ("symmetric", "zero_mean"),
            ("nonpositive", "less_power_only"),
        ]:
            result = solve_exact_map_noise_scipy(channel, 1, target_log2, kind, radius)
            rows.append({
                "secret_bits": bits,
                "target_log2": target_log2,
                "constraint": label,
                "status": result.status,
                "second_moment": result.objective,
                "rms_power_bins": math.sqrt(max(result.objective, 0.0)),
                "mean_noise": result.mean,
                "risk_log2": result.full_log_risk / math.log(2.0),
                "no_noise_log2": no_noise_log2,
                "prior_log2": prior_log2,
                "support_nonzero": int((result.probabilities > 1e-8).sum()),
            })

    write_csv(OUT / "power_focus_results.csv", rows)

    noise_series: Dict[str, List[Tuple[float, float]]] = {
        "zero-mean": [],
        "less-power only": [],
    }
    rms_series: Dict[str, List[Tuple[float, float]]] = {
        "zero-mean": [],
        "less-power only": [],
    }
    risk_series: Dict[str, List[Tuple[float, float]]] = {
        "no-noise exponent": [],
        "target exponent": [],
        "after optimization": [],
        "prior exponent": [],
    }

    for row in rows:
        label = "zero-mean" if row["constraint"] == "zero_mean" else "less-power only"
        noise_series[label].append((float(row["secret_bits"]), max(float(row["second_moment"]), 1e-6)))
        rms_series[label].append((float(row["secret_bits"]), max(float(row["rms_power_bins"]), 1e-6)))

    seen = set()
    for row in rows:
        bits = int(row["secret_bits"])
        if bits in seen:
            continue
        seen.add(bits)
        risk_series["no-noise exponent"].append((bits, -float(row["no_noise_log2"])))
        risk_series["target exponent"].append((bits, -target_log2))
        risk_series["prior exponent"].append((bits, -float(row["prior_log2"])))
    for row in rows:
        if row["constraint"] == "zero_mean":
            risk_series["after optimization"].append((float(row["secret_bits"]), -float(row["risk_log2"])))

    line_plot(
        OUT / "power_noise_second_moment_target247.png",
        "Aggregate AES Power: Noise for Target 2^-247",
        "Scalar Hamming-weight leakage; exact finite-sequence LP; zeros shown at floor",
        "Secret length l (bits)",
        "Noise second moment (power-bin^2)",
        {k: sorted(v) for k, v in noise_series.items()},
        log_y=True,
        y_floor=1e-6,
        x_ticks_override=lengths,
    )

    line_plot(
        OUT / "power_noise_rms_target247.png",
        "Aggregate AES Power: RMS Noise for Target 2^-247",
        "One power bin is one aggregate Hamming-weight count",
        "Secret length l (bits)",
        "RMS noise (power bins)",
        {k: sorted(v) for k, v in rms_series.items()},
        log_y=True,
        y_floor=1e-6,
        x_ticks_override=lengths,
    )

    line_plot(
        OUT / "power_posterior_exponents_target247.png",
        "Aggregate AES Power: Posterior Success Exponents",
        "Higher exponent means lower posterior success",
        "Secret length l (bits)",
        "-log2 posterior success",
        {k: sorted(v) for k, v in risk_series.items()},
        log_y=False,
        x_ticks_override=lengths,
        y_min_override=238.0,
        y_max_override=258.0,
    )

    print(f"wrote focused power results to {OUT}")


if __name__ == "__main__":
    run()
