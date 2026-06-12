"""Figures 5--8: entropy, mechanisms, and composition."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw

from convex_replot_experiments import font, line_plot
from discrete_noise_convex import (
    make_binomial_aggregate_channel,
    solve_exact_map_noise_scipy,
)
from prw_finite_alpha_power.run_power_prw_finite_alpha import (
    composition_prw_sigma2,
    uniform_prior,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "figures_5_to_8"
COMPOSITION_SEED = 20260610
COMPOSITION_TRIALS = 5
COMPOSITION_PRW_ALPHA = 256.0
COMPOSITION_ROUNDS = list(range(1, 51))
COMPOSITION_X_TICKS = [1, 10, 20, 30, 40, 50]
DISCRETE_NOISE_RADIUS = 2048
BOUNDED_LABEL = "Optimal zero-mean"
ONESIDED_LABEL = "Optimal non-negative noise"


def normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def gaussian_max_integral_uniform(bits: int, sigma: float) -> float:
    """Integral of max_h phi_sigma(y-h) for h=0,...,bits."""
    if sigma <= 0.0:
        return float(bits + 1)
    a = 0.5 / sigma
    phi = normal_cdf(a)
    return 2.0 * phi + (bits - 1.0) * (2.0 * phi - 1.0)


def gaussian_log2_success_uniform(bits: int, sigma: float) -> float:
    return -float(bits) + math.log2(gaussian_max_integral_uniform(bits, sigma))


def gaussian_sigma_for_target(bits: int, target_exponent: float) -> float:
    """Smallest sigma such that MAP success <= 2^{-target_exponent}."""
    target_integral = 2.0 ** (bits - target_exponent)
    if target_integral >= bits + 1:
        return 0.0
    if target_integral <= 1.0:
        return float("inf")

    def integral(sigma: float) -> float:
        return gaussian_max_integral_uniform(bits, sigma)

    lo, hi = 0.0, 1.0
    while integral(hi) > target_integral:
        hi *= 2.0
        if hi > 1e9:
            return float("inf")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if integral(mid) <= target_integral:
            hi = mid
        else:
            lo = mid
    return hi


def expected_abs_gaussian(sigma: float) -> float:
    if not math.isfinite(sigma):
        return float("inf")
    return sigma * math.sqrt(2.0 / math.pi)


def expected_abs_discrete(result) -> float:
    return float(np.sum(np.abs(result.support.astype(float)) * result.probabilities))


def second_moment_discrete(result) -> float:
    return float(np.sum((result.support.astype(float) ** 2) * result.probabilities))


def has_feasible_distribution(result) -> bool:
    return result.status in {"optimal", "zero_noise_feasible"} and np.isfinite(result.objective)


def discrete_for_target(bits: int, target_exponent: float, kind: str, radius: int):
    channel = make_binomial_aggregate_channel(bits)
    return solve_exact_map_noise_scipy(
        channel,
        1,
        -float(target_exponent),
        kind,
        radius,
        objective="expected_abs",
    )


def interval_map_expected_abs(span: int, map_factor: float, kind: str, radius: int) -> Tuple[float, str]:
    """Exact MAP interval optimizer for aggregate HW leakage under uniform secrets.

    The MAP factor for additive noise over an interval leakage domain is
    1 + span / support_size.  Minimizing E|N| gives the centered interval for
    symmetric noise and the initial interval for one-sided noise.
    """
    if map_factor <= 1.0:
        return float("inf"), "infeasible_below_no_leakage_bound"
    no_noise_factor = float(span + 1)
    if map_factor >= no_noise_factor:
        return 0.0, "zero_noise_feasible"
    if kind == "symmetric":
        max_support_size = float(2 * radius + 1)
        min_factor = 1.0 + float(span) / max_support_size
        if map_factor < min_factor:
            return float("inf"), "infeasible_support_radius"
        return float(span) / (4.0 * (map_factor - 1.0)), "interval_map_closed_form"
    if kind == "onesided":
        max_support_size = float(radius + 1)
        min_factor = 1.0 + float(span) / max_support_size
        if map_factor < min_factor:
            return float("inf"), "infeasible_support_radius"
        return float(span) / (2.0 * (map_factor - 1.0)), "interval_map_closed_form"
    raise ValueError(kind)


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


Panel = Dict[str, object]


def figure5(rows: List[dict]) -> Panel:
    lengths = [253, 254, 255, 256]
    target_exponents = [float(x) for x in np.arange(248.0, 252.0 + 1e-9, 0.25)]
    series: Dict[str, List[Tuple[float, float]]] = {}
    for bits in lengths:
        pts = []
        for s in target_exponents:
            sigma = gaussian_sigma_for_target(bits, s)
            pts.append((float(s), sigma))
            rows.append({
                "figure": 5,
                "mechanism": "Gaussian",
                "secret_bits": bits,
                "target_exponent": s,
                "gaussian_sigma": sigma,
                "y_value": sigma,
                "quantity": "std",
            })
        series[f"l={bits}"] = pts
    line_plot(
        OUT / "fig5_entropy_gaussian_std.png",
        "Entropy Effect under Gaussian Randomization",
        "Aggregate Hamming-weight leakage; full identification",
        "-log2 target posterior success",
        "Gaussian noise std sigma",
        series,
        x_ticks_override=[248, 249, 250, 251, 252],
        log_y=False,
    )
    return {
        "title": "(e) Effect from various l-bit secret entropy",
        "xlabel": "",
        "ylabel": "Gaussian noise std σ",
        "series": series,
        "x_ticks": [248, 249, 250, 251, 252],
        "inline_labels": True,
    }


def figure6(rows: List[dict]) -> Panel:
    bits = 256
    target_exponents = [float(x) for x in np.arange(248.0, 254.0 + 1e-9, 0.25)]
    series: Dict[str, List[Tuple[float, float]]] = {
        "Gaussian": [],
        BOUNDED_LABEL: [],
        ONESIDED_LABEL: [],
    }
    for s in target_exponents:
        sigma = gaussian_sigma_for_target(bits, s)
        gauss_abs = expected_abs_gaussian(sigma)
        series["Gaussian"].append((float(s), gauss_abs))
        rows.append({
            "figure": 6,
            "mechanism": "Gaussian",
            "secret_bits": bits,
            "target_exponent": s,
            "gaussian_sigma": sigma,
            "expected_abs_noise": gauss_abs,
            "second_moment": sigma * sigma,
            "status": "closed_form",
        })

        map_factor = 2.0 ** (bits - s)
        sym_abs, sym_status = interval_map_expected_abs(bits, map_factor, "symmetric", DISCRETE_NOISE_RADIUS)
        series[BOUNDED_LABEL].append((float(s), sym_abs))
        rows.append({
            "figure": 6,
            "mechanism": "optimal_discrete_symmetric",
            "secret_bits": bits,
            "target_exponent": s,
            "expected_abs_noise": sym_abs,
            "objective_value": sym_abs,
            "optimization_objective": "expected_abs",
            "risk_log2": -float(s) if math.isfinite(sym_abs) else "",
            "status": sym_status,
            "noise_radius": DISCRETE_NOISE_RADIUS,
            "accounting": "map_alpha_infinity",
        })

        pos_abs, pos_status = interval_map_expected_abs(bits, map_factor, "onesided", DISCRETE_NOISE_RADIUS)
        series[ONESIDED_LABEL].append((float(s), pos_abs))
        rows.append({
            "figure": 6,
            "mechanism": "optimal_discrete_nonnegative",
            "secret_bits": bits,
            "target_exponent": s,
            "expected_abs_noise": pos_abs,
            "objective_value": pos_abs,
            "optimization_objective": "expected_abs",
            "risk_log2": -float(s) if math.isfinite(pos_abs) else "",
            "status": pos_status,
            "noise_radius": DISCRETE_NOISE_RADIUS,
            "accounting": "map_alpha_infinity",
        })

    line_plot(
        OUT / "fig6_randomization_mechanisms.png",
        "Randomization Mechanisms for a 256-Bit Secret",
        "Objective for discrete mechanisms is minimum second moment; y-axis reports E|N|",
        "-log2 target posterior success",
        "Expected absolute noise E|N|",
        series,
        x_ticks_override=[248, 249, 250, 251, 252, 253, 254],
        log_y=False,
    )
    return {
        "title": "(f) Randomization mechanisms",
        "xlabel": "",
        "ylabel": "Expected absolute perturbation",
        "series": series,
        "x_ticks": [248, 249, 250, 251, 252, 253, 254],
    }


def local_exponent_for_composition(secret_bits: int, total_target_exponent: float, rounds: int) -> float:
    # Exact-MAP multiplicative accounting with equal per-round factors:
    # local MAP target = 2^{-l} * 2^{(l-s)/T}.
    return secret_bits - (secret_bits - total_target_exponent) / rounds


def composition_figure(rows: List[dict], fig_no: int, total_target_exponent: float, path_name: str) -> Panel:
    bits = 256
    rounds = COMPOSITION_ROUNDS
    gaussian_prior = uniform_prior(bits)
    rng = np.random.default_rng(COMPOSITION_SEED + fig_no)
    series: Dict[str, List[Tuple[float, float]]] = {
        "Gaussian": [],
        BOUNDED_LABEL: [],
        ONESIDED_LABEL: [],
    }

    # We explicitly generate public XOR masks for every trial/round to match the
    # composition interpretation.  Under a uniform 256-bit secret, the marginal
    # law of HW(X xor v_t) is exactly Binomial(256, 1/2) for every fixed v_t, so
    # the optimized local mechanism is invariant to the sampled mask.
    query_hamming_weights = rng.integers(0, 2, size=(COMPOSITION_TRIALS, len(rounds), bits), dtype=np.uint8).sum(axis=2)

    for t in rounds:
        local_s = local_exponent_for_composition(bits, total_target_exponent, t)
        sigma2 = composition_prw_sigma2(
            gaussian_prior,
            COMPOSITION_PRW_ALPHA,
            -float(total_target_exponent),
            t,
            gaussian_prior.prior_success,
        )
        sigma = math.sqrt(sigma2) if math.isfinite(sigma2) else float("inf")
        sigma2 = sigma * sigma if math.isfinite(sigma) else float("inf")
        gauss_abs = expected_abs_gaussian(sigma)
        series["Gaussian"].append((float(t), gauss_abs))
        rows.append({
            "figure": fig_no,
            "mechanism": "Gaussian",
            "secret_bits": bits,
            "total_target_exponent": total_target_exponent,
            "rounds": t,
            "local_target_exponent": local_s,
            "gaussian_sigma": sigma,
            "expected_abs_noise_per_round": gauss_abs,
            "second_moment": sigma2,
            "status": "finite_alpha_prw_equal_split",
            "composition_alpha": COMPOSITION_PRW_ALPHA,
            "trials": COMPOSITION_TRIALS,
            "query_model": "HW(X xor V_t)",
            "reported_quantity": "prefix_average_expected_abs_noise",
            "avg_query_hamming_weight": float(np.mean(query_hamming_weights[:, :t])),
        })

        map_factor = 2.0 ** (bits - local_s)
        for label, mechanism, kind in [
            (BOUNDED_LABEL, "optimal_discrete_symmetric", "symmetric"),
            (ONESIDED_LABEL, "optimal_discrete_nonnegative", "onesided"),
        ]:
            value, status = interval_map_expected_abs(bits, map_factor, kind, DISCRETE_NOISE_RADIUS)
            series[label].append((float(t), value))
            rows.append({
                "figure": fig_no,
                "mechanism": mechanism,
                "secret_bits": bits,
                "total_target_exponent": total_target_exponent,
                "rounds": t,
                "local_target_exponent": local_s,
                "expected_abs_noise_per_round": value,
                "objective_value": value,
                "optimization_objective": "expected_abs",
                "status": status,
                "noise_radius": DISCRETE_NOISE_RADIUS,
                "accounting": "map_alpha_infinity_equal_split",
                "query_model": "HW(X xor V_t)",
                "reported_quantity": "prefix_average_expected_abs_noise",
                "avg_query_hamming_weight": float(np.mean(query_hamming_weights[:, :t])),
            })

    plot_series = {
        name: [(x, y) for x, y in pts if math.isfinite(y)]
        for name, pts in series.items()
    }
    line_plot(
        OUT / path_name,
        f"Composition with Total Target 2^-{int(total_target_exponent)}",
        f"Gaussian: finite-alpha PRW alpha={int(COMPOSITION_PRW_ALPHA)}; optimized discrete: MAP interval radius={DISCRETE_NOISE_RADIUS}",
        "Number of iterations",
        "Prefix-average expected absolute noise (1/t) sum_i E|N_i|",
        plot_series,
        x_ticks_override=COMPOSITION_X_TICKS,
        log_y=False,
    )
    return {
        "title": f"({chr(ord('a') + fig_no - 1)}) Composition, full identification of 256-bit by 2^-{int(total_target_exponent)}",
        "xlabel": "",
        "ylabel": "Average perturbation up to iteration t" if fig_no == 7 else "",
        "series": series,
        "x_ticks": COMPOSITION_X_TICKS,
    }


def panel_colors(series_names: Sequence[str]) -> Dict[str, Tuple[int, int, int]]:
    entropy_palette = [
        (0, 114, 178),
        (213, 94, 0),
        (0, 158, 115),
        (204, 121, 167),
    ]
    mechanism_palette = {
        "Gaussian": (86, 180, 233),
        BOUNDED_LABEL: (117, 112, 179),
        ONESIDED_LABEL: (166, 86, 40),
    }
    colors: Dict[str, Tuple[int, int, int]] = {}
    for idx, name in enumerate(series_names):
        colors[name] = mechanism_palette.get(name, entropy_palette[idx % len(entropy_palette)])
    return colors


def draw_horizontal_label(
    d: ImageDraw.ImageDraw,
    xy: Tuple[float, float],
    label: str,
    color: Tuple[int, int, int],
    text_font,
    *,
    scale: int,
) -> float:
    x, y = xy
    d.line([(x, y + 15 * scale), (x + 42 * scale, y + 15 * scale)], fill=color, width=5 * scale)
    d.text((x + 52 * scale, y), label, fill=(28, 28, 28), font=text_font)
    return x + 52 * scale + d.textlength(label, font=text_font) + 28 * scale


def draw_rotated_label(
    img: Image.Image,
    text: str,
    center_x: float,
    center_y: float,
    text_font,
    *,
    scale: int,
) -> None:
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    text_w = int(measure.textlength(text, font=text_font) + 24 * scale)
    tmp = Image.new("RGBA", (text_w, 70 * scale), (255, 255, 255, 0))
    td = ImageDraw.Draw(tmp)
    td.text((0, 0), text, fill=(25, 25, 25), font=text_font)
    bbox = tmp.getbbox()
    if bbox:
        tmp = tmp.crop(bbox)
    rot = tmp.rotate(90, expand=True)
    img.paste(rot, (int(center_x - rot.width / 2), int(center_y - rot.height / 2)), rot)


def plot_combined_panels(path: Path, panels: Sequence[Panel]) -> None:
    final_width, final_height = 5600, 1300
    scale = 2
    width, height = final_width * scale, final_height * scale
    margin = {
        "left": 200 * scale,
        "right": 60 * scale,
        "top": 215 * scale,
        "bottom": 165 * scale,
    }
    gap_x = 150 * scale
    panel_w = (width - margin["left"] - margin["right"] - 3 * gap_x) / 4.0
    panel_h = height - margin["top"] - margin["bottom"]

    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    title_font = font(40 * scale, True)
    tick_font = font(33 * scale)
    axis_font = font(50 * scale, True)
    legend_font = font(38 * scale)

    def draw_log2_xlabel(center_x: float, y: float) -> None:
        prefix = "log"
        suffix = " of target posterior success rate"
        sub_font = font(39 * scale, True)
        prefix_w = d.textlength(prefix, font=axis_font)
        sub_w = d.textlength("2", font=sub_font)
        suffix_w = d.textlength(suffix, font=axis_font)
        total_w = prefix_w + sub_w + suffix_w + 14 * scale
        x = center_x - total_w / 2
        d.text((x, y), prefix, fill=(25, 25, 25), font=axis_font)
        sub_x = x + prefix_w + 3 * scale
        d.text((sub_x, y + 26 * scale), "2", fill=(25, 25, 25), font=sub_font)
        d.text((sub_x + sub_w + 8 * scale, y), suffix, fill=(25, 25, 25), font=axis_font)

    mechanism_names = list(panels[1]["series"].keys())  # type: ignore[index, union-attr]
    mechanism_colors = panel_colors(mechanism_names)
    y = 42 * scale

    f_x0 = margin["left"] + 1 * (panel_w + gap_x)
    h_x1 = margin["left"] + 3 * (panel_w + gap_x) + panel_w
    mech_width = sum(
        52 * scale + d.textlength(name, font=legend_font) + 28 * scale
        for name in mechanism_names
    )
    x = f_x0 + (h_x1 - f_x0 - mech_width) / 2.0
    for name in mechanism_names:
        x = draw_horizontal_label(d, (x, y), name, mechanism_colors[name], legend_font, scale=scale)

    for panel_idx, panel in enumerate(panels):
        title = str(panel["title"])
        xlabel = str(panel["xlabel"])
        ylabel = str(panel["ylabel"])
        inline_labels = bool(panel.get("inline_labels", False))
        series: Dict[str, Sequence[Tuple[float, float]]] = panel["series"]  # type: ignore[assignment]
        x_ticks = list(panel["x_ticks"])  # type: ignore[arg-type]
        colors = panel_colors(list(series.keys()))

        x0 = margin["left"] + panel_idx * (panel_w + gap_x)
        y0 = margin["top"]
        x1 = x0 + panel_w
        y1 = y0 + panel_h
        all_x = [x for pts in series.values() for x, _ in pts]
        all_y = [y for pts in series.values() for _, y in pts if y is not None and math.isfinite(float(y))]
        xmin, xmax = min(all_x), max(all_x)
        ymin = min(0.0, min(all_y))
        ymax = max(all_y)
        if ymax <= ymin:
            ymax = ymin + 1.0
        ymax = ymin + 1.12 * (ymax - ymin)

        def xm(value: float) -> float:
            return x0 + (value - xmin) / (xmax - xmin) * panel_w

        def ym(value: float) -> float:
            return y0 + (ymax - value) / (ymax - ymin) * panel_h

        d.text((x0, y0 - 58 * scale), title, fill=(25, 25, 25), font=title_font)
        d.rectangle([x0, y0, x1, y1], fill=(250, 250, 250), outline=(178, 178, 178))

        for tick in x_ticks:
            px = xm(float(tick))
            d.line([(px, y0), (px, y1)], fill=(226, 226, 226))
            label = f"{tick:g}"
            tw = d.textlength(label, font=tick_font)
            d.text((px - tw / 2, y1 + 12 * scale), label, fill=(45, 45, 45), font=tick_font)

        y_ticks = np.linspace(ymin, ymax, 6)
        for tick in y_ticks:
            py = ym(float(tick))
            d.line([(x0, py), (x1, py)], fill=(226, 226, 226))
            if abs(tick) >= 100:
                label = f"{tick:.0f}"
            elif abs(tick) >= 10:
                label = f"{tick:.1f}".rstrip("0").rstrip(".")
            else:
                label = f"{tick:.2g}"
            tw = d.textlength(label, font=tick_font)
            d.text((x0 - tw - 11 * scale, py - 13 * scale), label, fill=(45, 45, 45), font=tick_font)

        d.line([(x0, y1), (x1, y1)], fill=(30, 30, 30), width=3 * scale)
        d.line([(x0, y0), (x0, y1)], fill=(30, 30, 30), width=3 * scale)

        for name, pts in series.items():
            clean_pts = [(float(x), float(y)) for x, y in pts if y is not None and math.isfinite(float(y))]
            xy = [(xm(x), ym(y)) for x, y in clean_pts]
            color = colors[name]
            for a, b in zip(xy, xy[1:]):
                d.line([a, b], fill=color, width=5 * scale)
            for px, py in xy:
                d.ellipse([px - 4 * scale, py - 4 * scale, px + 4 * scale, py + 4 * scale], fill=color, outline="white", width=1 * scale)

        if inline_labels:
            label_x = x0 + 22 * scale
            label_y = y0 + 20 * scale
            for idx, name in enumerate(series.keys()):
                color = colors[name]
                yy = label_y + idx * 43 * scale
                d.line(
                    [(label_x, yy + 16 * scale), (label_x + 42 * scale, yy + 16 * scale)],
                    fill=color,
                    width=5 * scale,
                )
                d.text((label_x + 52 * scale, yy), name, fill=(28, 28, 28), font=legend_font)

        if xlabel:
            tw = d.textlength(xlabel, font=axis_font)
            d.text((x0 + panel_w / 2 - tw / 2, height - 66 * scale), xlabel, fill=(25, 25, 25), font=axis_font)
        if ylabel:
            draw_rotated_label(img, ylabel, x0 - 120 * scale, y0 + panel_h / 2, axis_font, scale=scale)

    ef_center = margin["left"] + (2 * panel_w + gap_x) / 2.0
    gh_left = margin["left"] + 2 * (panel_w + gap_x)
    gh_center = gh_left + (2 * panel_w + gap_x) / 2.0
    draw_log2_xlabel(ef_center, height - 90 * scale)
    iter_label = "Composition iteration"
    tw = d.textlength(iter_label, font=axis_font)
    d.text((gh_center - tw / 2, height - 90 * scale), iter_label, fill=(25, 25, 25), font=axis_font)

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    img = img.resize((final_width, final_height), resampling)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    panels = []
    panels.append(figure5(rows))
    print("wrote figure 5", flush=True)
    panels.append(figure6(rows))
    print("wrote figure 6", flush=True)
    panels.append(composition_figure(rows, 7, 240.0, "fig7_composition_target240.png"))
    print("wrote figure 7", flush=True)
    panels.append(composition_figure(rows, 8, 220.0, "fig8_composition_target220.png"))
    print("wrote figure 8", flush=True)
    plot_combined_panels(OUT / "fig5_8_combined_1x4_v5.png", panels)
    print("wrote fig5_8_combined_1x4_v5.png", flush=True)
    write_csv(OUT / "figures_5_to_8_results.csv", rows)
    print(f"wrote results to {OUT}")


if __name__ == "__main__":
    run()
