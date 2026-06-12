"""Compare history-dependent and mechanism-agnostic composition accounting.

The figure uses Gaussian perturbations.  The history-dependent curve uses the
finite-alpha PRW local condition with an equal split of the PRW composition
budget.  The mechanism-agnostic curve uses the balanced Holder allocation
p_t = T from the alpha-divergence agnostic theorem, so the per-round marginal
constraint is evaluated at order alpha*T.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw
from scipy.special import gammaln, logsumexp

from convex_replot_experiments import font
from rsa_timing_empirical_figures import (
    SEED as RSA_SEED,
    TimingCase,
    make_fast_composition_query,
    make_fast_composition_case,
    make_timing_case,
)


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "history_vs_agnostic"
DESKTOP = Path("/Users/hsxiao/Desktop")

BITS = 256
ALPHA = 256.0
TARGET_EXPONENTS = (240.0, 220.0)
MAX_ROUNDS = 100
ROUND_VALUES = tuple(range(1, 11)) + tuple(range(20, MAX_ROUNDS + 1, 10))
TIMING_SAMPLE_COUNT = 100_000
TIMING_TRIALS = 5


def expected_abs_gaussian(sigma: float) -> float:
    if not math.isfinite(sigma):
        return float("inf")
    return sigma * math.sqrt(2.0 / math.pi)


def log_binom_counts(bits: int) -> np.ndarray:
    k = np.arange(bits + 1, dtype=float)
    return gammaln(bits + 1.0) - gammaln(k + 1.0) - gammaln(bits - k + 1.0)


def power_marginal(bits: int) -> Tuple[np.ndarray, np.ndarray]:
    values = np.arange(bits + 1, dtype=float)
    log_probs = log_binom_counts(bits) - bits * math.log(2.0)
    return values, np.exp(log_probs)


def grouped_marginal(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    atoms, counts = np.unique(np.rint(values).astype(int), return_counts=True)
    probs = counts.astype(float) / float(np.sum(counts))
    return atoms.astype(float), probs


def integration_grid(values: np.ndarray, sigma: float) -> Tuple[np.ndarray, float]:
    lo = float(np.min(values)) - 8.0 * sigma
    hi = float(np.max(values)) + 8.0 * sigma
    dx = max(0.05, min(0.5, 0.025 * sigma))
    count = int(math.ceil((hi - lo) / dx)) + 1
    y = np.linspace(lo, hi, count)
    return y, float((hi - lo) / (count - 1))


def gaussian_logpdf_grid(y: np.ndarray, means: np.ndarray, sigma: float) -> np.ndarray:
    return -0.5 * math.log(2.0 * math.pi * sigma * sigma) - (
        (y[:, None] - means[None, :]) ** 2
    ) / (2.0 * sigma * sigma)


def log_optimized_prw_divergence(
    values: np.ndarray,
    probs: np.ndarray,
    bits: int,
    alpha_weight: float,
    divergence_order: float,
    sigma: float,
) -> float:
    """Optimized-reference log E[pi(X)^(alpha_weight-1) D_order(P_x||W)].

    For grouped uniform secrets, probs is the marginal mass of each leakage
    atom.  The optimal reference has density proportional to
    (sum_x pi(x)^alpha_weight p_x(y)^order)^(1/order).
    """
    if sigma <= 0.0:
        return float("inf")
    probs = np.asarray(probs, dtype=float)
    values = np.asarray(values, dtype=float)
    positive = probs > 0.0
    values = values[positive]
    probs = probs[positive]
    y, dx = integration_grid(values, sigma)
    logpdf = gaussian_logpdf_grid(y, values, sigma)
    log_prefix = -(alpha_weight - 1.0) * float(bits) * math.log(2.0)
    log_a = log_prefix + logsumexp(
        np.log(probs)[None, :] + divergence_order * logpdf,
        axis=1,
    )
    log_integral = float(logsumexp(log_a / divergence_order) + math.log(dx))
    return float(divergence_order * log_integral)


def solve_sigma_for_log_target(
    values: np.ndarray,
    probs: np.ndarray,
    bits: int,
    alpha_weight: float,
    divergence_order: float,
    target_log_value: float,
) -> float:
    log_no_leakage = -(alpha_weight - 1.0) * float(bits) * math.log(2.0)
    if log_no_leakage > target_log_value:
        return float("inf")

    def log_value(sigma: float) -> float:
        return log_optimized_prw_divergence(
            values,
            probs,
            bits,
            alpha_weight,
            divergence_order,
            sigma,
        )

    if log_value(0.25) <= target_log_value:
        return 0.0

    lo, hi = 0.25, 1.0
    while log_value(hi) > target_log_value:
        hi *= 2.0
        if hi > 1e9:
            return float("inf")

    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if log_value(mid) <= target_log_value:
            hi = mid
        else:
            lo = mid
    return hi


def history_target_log(bits: int, alpha: float, target_exponent: float, rounds: int) -> float:
    log_base = -(alpha - 1.0) * float(bits) * math.log(2.0)
    total_log_target = alpha * (-float(target_exponent)) * math.log(2.0)
    return log_base + (total_log_target - log_base) / float(rounds)


def agnostic_target_log(alpha: float, target_exponent: float) -> float:
    return alpha * (-float(target_exponent)) * math.log(2.0)


def curves_for_marginals(
    marginals_by_round: Sequence[Tuple[np.ndarray, np.ndarray]],
    bits: int,
    alpha: float,
    target_exponent: float,
    rounds: int,
) -> Tuple[float, float]:
    hist_target = history_target_log(bits, alpha, target_exponent, rounds)
    agn_target = agnostic_target_log(alpha, target_exponent)
    hist_values: List[float] = []
    agn_values: List[float] = []
    for values, probs in marginals_by_round:
        hist_sigma = solve_sigma_for_log_target(
            values,
            probs,
            bits,
            alpha,
            alpha,
            hist_target,
        )
        agn_sigma = solve_sigma_for_log_target(
            values,
            probs,
            bits,
            alpha,
            alpha * float(rounds),
            agn_target,
        )
        hist_values.append(expected_abs_gaussian(hist_sigma))
        agn_values.append(expected_abs_gaussian(agn_sigma))
    return float(np.mean(hist_values)), float(np.mean(agn_values))


def make_power_curves(rows: List[dict]) -> Dict[str, List[Tuple[float, float]]]:
    values, probs = power_marginal(BITS)
    marginals = [(values, probs)]
    series: Dict[str, List[Tuple[float, float]]] = {}
    for target in TARGET_EXPONENTS:
        hist_name = f"History-dependent, 2^-{int(target)}"
        agn_name = f"Agnostic, 2^-{int(target)}"
        series[hist_name] = []
        series[agn_name] = []
        for rounds in ROUND_VALUES:
            hist_abs, agn_abs = curves_for_marginals(
                marginals,
                BITS,
                ALPHA,
                target,
                rounds,
            )
            series[hist_name].append((float(rounds), hist_abs))
            series[agn_name].append((float(rounds), agn_abs))
            rows.append({
                "case": "power",
                "target_exponent": target,
                "rounds": rounds,
                "alpha": ALPHA,
                "method": "history_dependent",
                "expected_abs_gaussian": hist_abs,
                "divergence_order": ALPHA,
                "samples": "closed_form_binomial",
            })
            rows.append({
                "case": "power",
                "target_exponent": target,
                "rounds": rounds,
                "alpha": ALPHA,
                "method": "mechanism_agnostic",
                "expected_abs_gaussian": agn_abs,
                "divergence_order": ALPHA * float(rounds),
                "samples": "closed_form_binomial",
            })
    return series


def timing_case_marginal(case: TimingCase) -> Tuple[np.ndarray, np.ndarray]:
    probs = case.leakage_counts.astype(float) / float(np.sum(case.leakage_counts))
    return case.leakage_values.astype(float), probs


def precompute_timing_representative_marginals() -> List[Tuple[np.ndarray, np.ndarray]]:
    base256 = make_timing_case(BITS, label="Uniform 256-bit", seed=RSA_SEED + 2)
    comp_rng = np.random.default_rng(RSA_SEED + 1000)
    marginals: List[Tuple[np.ndarray, np.ndarray]] = []
    for trial_idx in range(TIMING_TRIALS):
        query = make_fast_composition_query(comp_rng, BITS)
        case = make_fast_composition_case(
            base256.samples,
            bits=BITS,
            label=f"representative timing query {trial_idx + 1}",
            query=query,
            log_prior=base256.log_prior,
            prior_success_log2=base256.prior_success_log2,
        )
        marginals.append(timing_case_marginal(case))
    return marginals


def make_timing_curves(rows: List[dict]) -> Dict[str, List[Tuple[float, float]]]:
    timing_marginals = precompute_timing_representative_marginals()
    series: Dict[str, List[Tuple[float, float]]] = {}
    for target in TARGET_EXPONENTS:
        hist_name = f"History-dependent, 2^-{int(target)}"
        agn_name = f"Agnostic, 2^-{int(target)}"
        series[hist_name] = []
        series[agn_name] = []
        for rounds in ROUND_VALUES:
            hist_trial_values = []
            agn_trial_values = []
            for trial_idx, marginal in enumerate(timing_marginals, start=1):
                hist_abs, agn_abs = curves_for_marginals(
                    [marginal],
                    BITS,
                    ALPHA,
                    target,
                    rounds,
                )
                hist_trial_values.append(hist_abs)
                agn_trial_values.append(agn_abs)
                rows.append({
                    "case": "timing",
                    "target_exponent": target,
                    "rounds": rounds,
                    "trial": trial_idx,
                    "alpha": ALPHA,
                    "method": "history_dependent",
                    "expected_abs_gaussian": hist_abs,
                    "divergence_order": ALPHA,
                    "samples": TIMING_SAMPLE_COUNT,
                    "query_model": "representative_random_query_marginal",
                })
                rows.append({
                    "case": "timing",
                    "target_exponent": target,
                    "rounds": rounds,
                    "trial": trial_idx,
                    "alpha": ALPHA,
                    "method": "mechanism_agnostic",
                    "expected_abs_gaussian": agn_abs,
                    "divergence_order": ALPHA * float(rounds),
                    "samples": TIMING_SAMPLE_COUNT,
                    "query_model": "representative_random_query_marginal",
                })
            series[hist_name].append((float(rounds), float(np.mean(hist_trial_values))))
            series[agn_name].append((float(rounds), float(np.mean(agn_trial_values))))
    return series


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


def draw_figure(
    power_series: Dict[str, List[Tuple[float, float]]],
    timing_series: Dict[str, List[Tuple[float, float]]],
) -> None:
    palettes = {
        "History-based": "#0072B2",
        "Mechanism-agnostic": "#D55E00",
    }
    panels: List[Tuple[str, Dict[str, List[Tuple[float, float]]]]] = []
    for prefix, source, letter_start in [
        ("Power leakage", power_series, "a"),
        ("RSA timing leakage", timing_series, "c"),
    ]:
        for offset, target in enumerate(TARGET_EXPONENTS):
            letter = chr(ord(letter_start) + offset)
            hist_name = f"History-dependent, 2^-{int(target)}"
            agn_name = f"Agnostic, 2^-{int(target)}"
            panels.append((
                f"({letter}) {prefix}, target 2^-{int(target)}",
                {
                    "History-based": source[hist_name],
                    "Mechanism-agnostic": source[agn_name],
                },
            ))
    OUT.mkdir(parents=True, exist_ok=True)
    combined = render_panels(panels, palettes, combined=True)
    combined.save(OUT / "history_vs_agnostic_composition.png")
    combined.save(DESKTOP / "history_vs_agnostic_composition.png")

    for name, series in [
        ("power_history_vs_agnostic_composition.png", {
            "History-based": power_series[f"History-dependent, 2^-{int(TARGET_EXPONENTS[0])}"],
            "Mechanism-agnostic": power_series[f"Agnostic, 2^-{int(TARGET_EXPONENTS[0])}"],
        }),
        ("timing_history_vs_agnostic_composition.png", {
            "History-based": timing_series[f"History-dependent, 2^-{int(TARGET_EXPONENTS[0])}"],
            "Mechanism-agnostic": timing_series[f"Agnostic, 2^-{int(TARGET_EXPONENTS[0])}"],
        }),
    ]:
        img = render_panels([(name.replace("_", " ").replace(".png", ""), series)], palettes, combined=False)
        img.save(OUT / name)
        img.save(DESKTOP / name)


def hex_to_rgb(value: str) -> Tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))  # type: ignore[return-value]


def draw_dashed_line(
    d: ImageDraw.ImageDraw,
    points: Sequence[Tuple[float, float]],
    *,
    fill: Tuple[int, int, int],
    width: int,
    dash: int,
    gap: int,
) -> None:
    for start, end in zip(points, points[1:]):
        x0, y0 = start
        x1, y1 = end
        dx = x1 - x0
        dy = y1 - y0
        dist = math.hypot(dx, dy)
        if dist == 0:
            continue
        ux, uy = dx / dist, dy / dist
        pos = 0.0
        while pos < dist:
            seg_end = min(dist, pos + dash)
            d.line(
                [
                    (x0 + ux * pos, y0 + uy * pos),
                    (x0 + ux * seg_end, y0 + uy * seg_end),
                ],
                fill=fill,
                width=width,
            )
            pos += dash + gap


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


def format_tick(value: float) -> str:
    if value >= 1000:
        return f"{value:.0f}"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2g}"


def render_panels(
    panels: Sequence[Tuple[str, Dict[str, List[Tuple[float, float]]]]],
    palettes: Dict[str, str],
    *,
    combined: bool,
) -> Image.Image:
    scale = 2
    final_width = 6000 if combined else 1900
    final_height = 1420 if combined else 1220
    width, height = final_width * scale, final_height * scale
    margin = {
        "left": 260 * scale,
        "right": 65 * scale,
        "top": 305 * scale if combined else 200 * scale,
        "bottom": 185 * scale,
    }
    gap_x = 105 * scale if combined else 150 * scale
    panel_count = len(panels)
    panel_w = (width - margin["left"] - margin["right"] - gap_x * (panel_count - 1)) / panel_count
    panel_h = height - margin["top"] - margin["bottom"]
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    title_font = font(56 * scale, True)
    axis_font = font(56 * scale, True)
    tick_font = font(56 * scale)
    legend_font = font(56 * scale)
    small_font = font(56 * scale)
    x_ticks = [1, 20, 40, 60, 80, MAX_ROUNDS]
    marker_every = max(1, MAX_ROUNDS // 20)

    if combined:
        names = list(panels[0][1].keys())
        legend_y = 38 * scale
        total_w = sum(78 * scale + d.textlength(name.replace("$", ""), font=legend_font) + 48 * scale for name in names)
        x = (width - total_w) / 2
        for name in names:
            color = hex_to_rgb(palettes[name])
            y_line = legend_y + 24 * scale
            if name.startswith("Mechanism"):
                draw_dashed_line(d, [(x, y_line), (x + 66 * scale, y_line)], fill=color, width=7 * scale, dash=16 * scale, gap=9 * scale)
            else:
                d.line([(x, y_line), (x + 66 * scale, y_line)], fill=color, width=7 * scale)
            text = name
            d.text((x + 82 * scale, legend_y), text, fill=(25, 25, 25), font=legend_font)
            x += 78 * scale + d.textlength(text, font=legend_font) + 48 * scale

    for panel_idx, (title, series) in enumerate(panels):
        x0 = margin["left"] + panel_idx * (panel_w + gap_x)
        y0 = margin["top"]
        x1 = x0 + panel_w
        y1 = y0 + panel_h
        all_y = [y for pts in series.values() for _, y in pts if math.isfinite(y)]
        ymin = 0.0
        ymax = max(all_y) if all_y else 1.0
        ymax = ymax * 1.10 if ymax > 0 else 1.0

        def xm(value: float) -> float:
            return x0 + (value - 1.0) / (MAX_ROUNDS - 1.0) * panel_w

        def ym(value: float) -> float:
            return y0 + (ymax - value) / (ymax - ymin) * panel_h

        d.text((x0, y0 - 64 * scale), title, fill=(25, 25, 25), font=title_font)
        d.rectangle([x0, y0, x1, y1], fill=(250, 250, 250), outline=(178, 178, 178))
        for tick in x_ticks:
            px = xm(float(tick))
            d.line([(px, y0), (px, y1)], fill=(226, 226, 226))
            label = str(tick)
            tw = d.textlength(label, font=tick_font)
            d.text((px - tw / 2, y1 + 12 * scale), label, fill=(45, 45, 45), font=tick_font)
        for tick in np.linspace(ymin, ymax, 6):
            py = ym(float(tick))
            d.line([(x0, py), (x1, py)], fill=(226, 226, 226))
            label = format_tick(float(tick))
            tw = d.textlength(label, font=tick_font)
            d.text((x0 - tw - 10 * scale, py - 13 * scale), label, fill=(45, 45, 45), font=tick_font)
        d.line([(x0, y1), (x1, y1)], fill=(30, 30, 30), width=3 * scale)
        d.line([(x0, y0), (x0, y1)], fill=(30, 30, 30), width=3 * scale)

        for name, pts in series.items():
            color = hex_to_rgb(palettes[name])
            mapped = [(xm(x), ym(y)) for x, y in pts]
            if name.startswith("Mechanism"):
                draw_dashed_line(d, mapped, fill=color, width=5 * scale, dash=13 * scale, gap=8 * scale)
            else:
                for a, b in zip(mapped, mapped[1:]):
                    d.line([a, b], fill=color, width=5 * scale)
            for idx, (px, py) in enumerate(mapped):
                if idx % marker_every != 0 and idx != len(mapped) - 1:
                    continue
                d.ellipse([px - 5 * scale, py - 5 * scale, px + 5 * scale, py + 5 * scale], fill=color, outline="white", width=1 * scale)

        xlabel = "Composition iteration"
        tw = d.textlength(xlabel, font=axis_font)
        d.text((x0 + panel_w / 2 - tw / 2, height - 92 * scale), xlabel, fill=(25, 25, 25), font=axis_font)
        if panel_idx == 0:
            draw_rotated_label(img, "Expected absolute Gaussian perturbation", 82 * scale, y0 + panel_h / 2, axis_font, scale=scale)
        if not combined:
            legend_x = x0 + 20 * scale
            legend_y = y0 + 20 * scale
            for idx, name in enumerate(series.keys()):
                color = hex_to_rgb(palettes[name])
                yy = legend_y + idx * 42 * scale
                if name.startswith("Mechanism"):
                    draw_dashed_line(d, [(legend_x, yy + 16 * scale), (legend_x + 48 * scale, yy + 16 * scale)], fill=color, width=5 * scale, dash=12 * scale, gap=7 * scale)
                else:
                    d.line([(legend_x, yy + 16 * scale), (legend_x + 48 * scale, yy + 16 * scale)], fill=color, width=5 * scale)
                d.text((legend_x + 60 * scale, yy), name, fill=(25, 25, 25), font=small_font)

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    return img.resize((final_width, final_height), resampling)


def run() -> None:
    rows: List[dict] = []
    power_series = make_power_curves(rows)
    print("computed power curves", flush=True)
    timing_series = make_timing_curves(rows)
    print("computed timing curves", flush=True)
    draw_figure(power_series, timing_series)
    write_csv(OUT / "history_vs_agnostic_composition.csv", rows)
    print(f"wrote {OUT / 'history_vs_agnostic_composition.png'}", flush=True)
    print(f"wrote {DESKTOP / 'history_vs_agnostic_composition.png'}", flush=True)


if __name__ == "__main__":
    run()
