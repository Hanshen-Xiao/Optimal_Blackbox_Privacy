"""Replot finite-domain experiments with convex sequence noise optimization.

This script uses the updated experiment setting:

* no confidence corrections;
* finite empirical/population distributions are used directly;
* AES power is the one-dimensional aggregate Hamming-weight channel;
* the noise distribution is optimized as a probability sequence on a finite
  support, minimizing E[N^2].

The old one-parameter discrete-Gaussian search is kept only as a baseline for
comparison.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from discrete_leakage_experiments import (
    calibrate_noise,
    make_aes_power_channel,
    make_rsa_timing_bit_channel,
)
from discrete_noise_convex import (
    make_binomial_aggregate_channel,
    solve_exact_map_noise,
    solve_prw_alpha_noise,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "convex"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    paths = []
    if bold:
        paths.extend([
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
            "/Library/Fonts/Arial Bold.ttf",
        ])
    paths.extend([
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
    ])
    for p in paths:
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()


def line_plot(
    path: Path,
    title: str,
    subtitle: str,
    xlabel: str,
    ylabel: str,
    series: Dict[str, Sequence[Tuple[float, float]]],
    *,
    log_x: bool = False,
    log_y: bool = False,
    y_floor: float = 1e-9,
    x_ticks_override: Sequence[float] | None = None,
    y_min_override: float | None = None,
    y_max_override: float | None = None,
) -> None:
    w, h = 1450, 900
    margin = {"left": 150, "right": 300, "top": 120, "bottom": 125}
    pw, ph = w - margin["left"] - margin["right"], h - margin["top"] - margin["bottom"]
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    ft, fs, fa, fk = font(32, True), font(18), font(23), font(18)
    d.text((margin["left"], 34), title, fill=(25, 25, 25), font=ft)
    d.text((margin["left"], 78), subtitle, fill=(80, 80, 80), font=fs)

    all_x = [x for pts in series.values() for x, _ in pts]
    all_y = [y for pts in series.values() for _, y in pts]
    xmin, xmax = min(all_x), max(all_x)
    if xmin == xmax:
        xmin -= 1.0
        xmax += 1.0
    if log_y:
        positive = [max(y, y_floor) for y in all_y]
        ymin = 10 ** math.floor(math.log10(max(min(positive), y_floor)))
        ymax = 10 ** math.ceil(math.log10(max(positive)))
        if ymin == ymax:
            ymax *= 10.0
    else:
        ymin = min(0.0, min(all_y))
        ymax = max(all_y)
        if ymax <= ymin:
            ymax = ymin + 1.0
        ymax = ymin + 1.12 * (ymax - ymin)
    if y_min_override is not None:
        ymin = y_min_override
    if y_max_override is not None:
        ymax = y_max_override

    def xm(x: float) -> float:
        if log_x:
            return margin["left"] + (math.log10(x) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin)) * pw
        return margin["left"] + (x - xmin) / (xmax - xmin) * pw

    def ym(y: float) -> float:
        if log_y:
            y = max(y, y_floor)
            return margin["top"] + (math.log10(ymax) - math.log10(y)) / (math.log10(ymax) - math.log10(ymin)) * ph
        return margin["top"] + (ymax - y) / (ymax - ymin) * ph

    d.rectangle(
        [margin["left"], margin["top"], margin["left"] + pw, margin["top"] + ph],
        fill=(250, 250, 250),
        outline=(180, 180, 180),
    )

    if x_ticks_override is not None:
        x_ticks = list(x_ticks_override)
    elif log_x:
        x_ticks = []
        for exp in range(int(math.floor(math.log10(xmin))), int(math.ceil(math.log10(xmax))) + 1):
            for mult in [1, 2, 5]:
                val = mult * 10**exp
                if xmin <= val <= xmax:
                    x_ticks.append(val)
    else:
        x_ticks = [xmin + (xmax - xmin) * i / 5 for i in range(6)]
    for x in x_ticks:
        px = xm(x)
        d.line([(px, margin["top"]), (px, margin["top"] + ph)], fill=(226, 226, 226))
        label = f"{x:g}"
        tw = d.textlength(label, font=fk)
        d.text((px - tw / 2, margin["top"] + ph + 14), label, fill=(45, 45, 45), font=fk)

    if log_y:
        y_ticks = []
        for exp in range(int(math.floor(math.log10(ymin))), int(math.ceil(math.log10(ymax))) + 1):
            for mult in [1, 2, 5]:
                val = mult * 10**exp
                if ymin <= val <= ymax:
                    y_ticks.append(val)
    else:
        y_ticks = [ymin + (ymax - ymin) * i / 5 for i in range(6)]
    for y in y_ticks:
        py = ym(y)
        d.line([(margin["left"], py), (margin["left"] + pw, py)], fill=(226, 226, 226))
        if log_y:
            label = f"{y:.2g}" if abs(y) < 1000 else f"{y/1000:.1f}k"
        elif abs(y) >= 1000:
            label = f"{y/1000:.1f}k"
        elif abs(y) >= 100:
            label = f"{y:.0f}"
        elif abs(y) >= 10:
            label = f"{y:.1f}".rstrip("0").rstrip(".")
        else:
            label = f"{y:.2g}"
        tw = d.textlength(label, font=fk)
        d.text((margin["left"] - tw - 12, py - 9), label, fill=(45, 45, 45), font=fk)

    d.line([(margin["left"], margin["top"] + ph), (margin["left"] + pw, margin["top"] + ph)], fill=(30, 30, 30), width=2)
    d.line([(margin["left"], margin["top"]), (margin["left"], margin["top"] + ph)], fill=(30, 30, 30), width=2)

    colors = [
        (0, 114, 178),
        (213, 94, 0),
        (0, 158, 115),
        (204, 121, 167),
        (86, 180, 233),
        (230, 159, 0),
        (90, 90, 90),
    ]
    for idx, (name, pts) in enumerate(series.items()):
        color = colors[idx % len(colors)]
        xy = [(xm(x), ym(y)) for x, y in pts]
        for a, b in zip(xy, xy[1:]):
            d.line([a, b], fill=color, width=4)
        for px, py in xy:
            d.ellipse([px - 6, py - 6, px + 6, py + 6], fill=color, outline="white", width=2)
        lx = margin["left"] + pw + 35
        ly = margin["top"] + 20 + idx * 34
        d.line([(lx, ly + 10), (lx + 45, ly + 10)], fill=color, width=4)
        d.text((lx + 58, ly), name, fill=(35, 35, 35), font=fk)

    tw = d.textlength(xlabel, font=fa)
    d.text((margin["left"] + pw / 2 - tw / 2, h - 72), xlabel, fill=(30, 30, 30), font=fa)
    label_img = Image.new("RGBA", (720, 40), (255, 255, 255, 0))
    ld = ImageDraw.Draw(label_img)
    ld.text((0, 0), ylabel, fill=(30, 30, 30), font=fa)
    rot = label_img.rotate(90, expand=True)
    img.paste(rot, (30, margin["top"] + ph // 2 - rot.height // 2), rot)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


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


def target_for_composition(secret_bits: int, alpha: float, total_target_log2: float, rounds: int) -> float:
    log_base = -secret_bits * (alpha - 1.0) * math.log(2.0)
    log_target_alpha = alpha * total_target_log2 * math.log(2.0)
    log_r_round = (log_target_alpha - log_base) / rounds
    return (log_base + log_r_round) / alpha / math.log(2.0)


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    entropies = [208, 224, 240, 256]
    targets = [-200, -150, -100]
    radius = 64

    rsa = make_rsa_timing_bit_channel(unit_ns=1.0)
    rsa_biased = make_rsa_timing_bit_channel(bit_prob_one=0.6, unit_ns=1.0)
    aes_byte = make_aes_power_channel()

    # 1. AES aggregate scalar leakage.
    for l in entropies:
        aes_agg = make_binomial_aggregate_channel(l)
        for target in targets:
            res = solve_exact_map_noise(aes_agg, 1, target, "symmetric", radius)
            rows.append({
                "experiment": "aes_aggregate",
                "channel": "AES aggregate HW",
                "entropy_bits": l,
                "target_log2": target,
                "constraint": "zero_mean",
                "method": "Exact sequence LP",
                "status": res.status,
                "second_moment": res.objective,
                "mean_noise": res.mean,
                "risk_log2": res.full_log_risk / math.log(2.0),
            })

    # 2. AES model comparison for the strongest requested target.
    for l in entropies:
        aes_agg = make_binomial_aggregate_channel(l)
        agg = solve_exact_map_noise(aes_agg, 1, -200, "symmetric", radius)
        byte_lp = solve_exact_map_noise(aes_byte, l // 8, -200, "symmetric", radius)
        byte_param = calibrate_noise(aes_byte, l // 8, -200, "symmetric", "prw", 128.0, 512, 600.0)
        for label, value, status in [
            ("Aggregate exact LP", agg.objective, agg.status),
            ("Byte-vector exact LP", byte_lp.objective, byte_lp.status),
            ("Byte-vector Gaussian family", byte_param.second_moment, byte_param.status),
        ]:
            rows.append({
                "experiment": "aes_model_comparison",
                "channel": label,
                "entropy_bits": l,
                "target_log2": -200,
                "constraint": "zero_mean",
                "method": label,
                "status": status,
                "second_moment": value,
            })

    # 3. RSA exact sequence LP versus the old Gaussian-family baseline.
    for target in targets:
        for kind, label in [("symmetric", "zero_mean"), ("onesided", "delay_only")]:
            exact = solve_exact_map_noise(rsa, 256, target, kind, radius)
            gaussian = calibrate_noise(rsa, 256, target, kind, "prw", 512.0, 512, 600.0)
            for method, objective, mean, status in [
                ("Exact sequence LP", exact.objective, exact.mean, exact.status),
                ("Gaussian-family PRW512", gaussian.second_moment, gaussian.mean, gaussian.status),
            ]:
                rows.append({
                    "experiment": "rsa_constraint",
                    "channel": "RSA timing",
                    "entropy_bits": 256,
                    "target_log2": target,
                    "constraint": label,
                    "method": method,
                    "status": status,
                    "second_moment": objective,
                    "rms_ns_per_bit": math.sqrt(max(objective, 0.0)),
                    "mean_delay_ns_per_bit": mean,
                    "mean_delay_ns_total": mean * 256,
                })

    # 4. Alpha convergence for convex PRW sequence optimization.
    exact_ref = solve_exact_map_noise(rsa, 256, -150, "symmetric", radius)
    for alpha in [2, 4, 8, 16, 32, 64, 128]:
        res = solve_prw_alpha_noise(rsa, 256, -150, "symmetric", 16, float(alpha))
        rows.append({
            "experiment": "alpha_convergence",
            "channel": "RSA timing",
            "entropy_bits": 256,
            "target_log2": -150,
            "constraint": "zero_mean",
            "method": f"Convex PRW alpha={alpha}",
            "alpha": alpha,
            "status": res.status,
            "second_moment": res.objective if res.status in {"optimal", "optimal_inaccurate", "zero_noise_feasible"} else "",
        })
    rows.append({
        "experiment": "alpha_convergence",
        "channel": "RSA timing",
        "entropy_bits": 256,
        "target_log2": -150,
        "constraint": "zero_mean",
        "method": "Exact MAP LP",
        "alpha": "exact",
        "status": exact_ref.status,
        "second_moment": exact_ref.objective,
    })

    # 5. PAC/KL comparison, using non-uniform RSA as the second case.
    comparison_cases = [
        (rsa, "uniform", -150),
        (rsa_biased, "nonuniform p=0.6", -120),
    ]
    for channel, prior, target in comparison_cases:
        exact = solve_exact_map_noise(channel, 256, target, "symmetric", radius)
        prw = solve_prw_alpha_noise(channel, 256, target, "symmetric", 16, 64.0)
        gaussian = calibrate_noise(channel, 256, target, "symmetric", "prw", 512.0, 512, 600.0)
        pac = calibrate_noise(channel, 256, target, "symmetric", "kl", None, 512, 600.0)
        for method, objective, status in [
            ("Exact sequence LP", exact.objective, exact.status),
            ("Convex PRW alpha=64", prw.objective, prw.status),
            ("Gaussian-family PRW512", gaussian.second_moment, gaussian.status),
            ("PAC KL/MI Gaussian", pac.second_moment if pac.status != "not_reached" else 1e5, pac.status),
        ]:
            rows.append({
                "experiment": "kl_comparison",
                "prior": prior,
                "channel": "RSA timing",
                "target_log2": target,
                "constraint": "zero_mean",
                "method": method,
                "status": status,
                "second_moment": objective,
            })

    # 6. Composition via equal PRW alpha budget split.
    alpha = 64.0
    total_target = -150.0
    for rounds in range(1, 9):
        local_target = target_for_composition(256, alpha, total_target, rounds)
        convex = solve_prw_alpha_noise(rsa, 256, local_target, "onesided", radius, alpha)
        gaussian = calibrate_noise(rsa, 256, local_target, "onesided", "prw", alpha, 512, 600.0)
        for method, objective, mean, status in [
            ("Convex PRW alpha=64", convex.objective, convex.mean, convex.status),
            ("Gaussian-family PRW64", gaussian.second_moment, gaussian.mean, gaussian.status),
        ]:
            rows.append({
                "experiment": "composition",
                "channel": "RSA timing",
                "entropy_bits": 256,
                "target_log2": total_target,
                "local_target_log2": local_target,
                "constraint": "delay_only",
                "method": method,
                "T": rounds,
                "status": status,
                "avg_second_moment_per_iter": objective,
                "mean_delay_ns_per_iter": mean * 256,
            })

    write_csv(OUT / "convex_experiment_results.csv", rows)

    aes_series: Dict[str, List[Tuple[float, float]]] = {}
    for target in targets:
        pts = [
            (float(r["entropy_bits"]), float(r["second_moment"]))
            for r in rows
            if r.get("experiment") == "aes_aggregate" and r.get("target_log2") == target
        ]
        aes_series[f"target 2^{target}"] = sorted(pts)
    line_plot(
        OUT / "aes_aggregate_entropy_targets.png",
        "AES Aggregate Power: Requested Targets Need No Added Noise",
        "Scalar aggregate Hamming-weight leakage, exact finite-sequence LP",
        "Secret entropy l (bits)",
        "Noise second moment",
        aes_series,
        log_y=False,
    )

    aes_cmp: Dict[str, List[Tuple[float, float]]] = {}
    for method in ["Aggregate exact LP", "Byte-vector exact LP", "Byte-vector Gaussian family"]:
        aes_cmp[method] = sorted([
            (float(r["entropy_bits"]), float(r["second_moment"]))
            for r in rows
            if r.get("experiment") == "aes_model_comparison" and r.get("method") == method
        ])
    line_plot(
        OUT / "aes_model_comparison.png",
        "AES Leakage Model Comparison",
        "Target 2^-200; aggregate scalar leakage is much weaker than byte-vector leakage",
        "Secret entropy l (bits)",
        "Noise second moment",
        aes_cmp,
        log_y=False,
    )

    rsa_cmp: Dict[str, List[Tuple[float, float]]] = {}
    for constraint, pretty in [("zero_mean", "zero-mean"), ("delay_only", "delay-only")]:
        for method, short in [("Exact sequence LP", "exact"), ("Gaussian-family PRW512", "Gaussian")]:
            rsa_cmp[f"{pretty} {short}"] = sorted([
                (-float(r["target_log2"]), max(float(r["second_moment"]), 1e-9))
                for r in rows
                if r.get("experiment") == "rsa_constraint"
                and r.get("constraint") == constraint
                and r.get("method") == method
            ])
    line_plot(
        OUT / "rsa_constraint_comparison.png",
        "RSA Timing: Convex Sequence Noise vs Gaussian Family",
        "Exact finite sequence LP substantially reduces E[N^2]",
        "Required posterior exponent s for target 2^-s",
        "Per-bit noise second moment (ns^2)",
        rsa_cmp,
        log_y=True,
    )

    alpha_pts = sorted([
        (float(r["alpha"]), max(float(r["second_moment"]), 1e-9))
        for r in rows
        if r.get("experiment") == "alpha_convergence"
        and isinstance(r.get("alpha"), int)
        and r.get("second_moment") != ""
    ])
    exact_pts = [(x, exact_ref.objective) for x, _ in alpha_pts]
    line_plot(
        OUT / "rsa_alpha_convergence_convex.png",
        "RSA Timing: Convex PRW Alpha Convergence",
        "Target 2^-150, zero-mean finite-sequence noise",
        "alpha",
        "Per-bit noise second moment (ns^2)",
        {"Convex PRW": alpha_pts, "Exact MAP LP": exact_pts},
        log_x=True,
        log_y=True,
    )

    method_index = {
        "Exact sequence LP": 1,
        "Convex PRW alpha=64": 2,
        "Gaussian-family PRW512": 3,
        "PAC KL/MI Gaussian": 4,
    }
    kl_series: Dict[str, List[Tuple[float, float]]] = {}
    for prior in ["uniform", "nonuniform p=0.6"]:
        kl_series[prior] = sorted([
            (float(method_index[str(r["method"])]), max(float(r["second_moment"]), 1e-9))
            for r in rows
            if r.get("experiment") == "kl_comparison" and r.get("prior") == prior
        ])
    line_plot(
        OUT / "rsa_pac_comparison_convex.png",
        "RSA Timing: PAC KL/MI Comparison",
        "Index: 1 Exact LP, 2 Convex PRW, 3 Gaussian PRW, 4 PAC KL/MI cap",
        "Method index",
        "Per-bit noise second moment (ns^2)",
        kl_series,
        log_y=True,
    )

    comp_series: Dict[str, List[Tuple[float, float]]] = {}
    for method in ["Convex PRW alpha=64", "Gaussian-family PRW64"]:
        comp_series[method] = sorted([
            (float(r["T"]), max(float(r["avg_second_moment_per_iter"]), 1e-9))
            for r in rows
            if r.get("experiment") == "composition" and r.get("method") == method
        ])
    line_plot(
        OUT / "rsa_composition_convex.png",
        "RSA Timing Composition",
        "Fixed total target 2^-150, one-sided delay noise, equal PRW alpha=64 split",
        "Composition rounds T",
        "Average per-round second moment (ns^2)",
        comp_series,
        log_y=True,
    )

    print(f"wrote convex experiment results to {OUT}")


if __name__ == "__main__":
    run()
