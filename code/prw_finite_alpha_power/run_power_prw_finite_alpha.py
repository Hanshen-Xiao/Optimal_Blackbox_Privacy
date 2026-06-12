"""Finite-alpha PRW Gaussian experiments for aggregate power leakage.

This script is intentionally self-contained.  It preserves the finite-alpha
Gaussian accounting profile used by the earlier AES notebook:

    alpha * log int (sum_h p_H(h) phi_sigma(y-h)^alpha)^{1/alpha} dy.

The current main experiment script uses an exact-MAP local target in the
composition panels.  This script instead splits the finite-alpha Bernoulli
divergence budget evenly across composed rounds and solves directly in sigma^2.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.optimize import brentq
from scipy.special import gammaln, logsumexp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
MAIN_ALPHA = 256.0
CURVE_ALPHAS = [16.0, 64.0, 256.0]
BUDGET_TABLES: Dict[Tuple[str, float], Tuple[np.ndarray, np.ndarray]] = {}


@dataclass
class AggregatePrior:
    label: str
    bits: int
    h: np.ndarray
    log_count: np.ndarray
    log_p_h: np.ndarray
    log_atom: np.ndarray

    @property
    def prior_success(self) -> float:
        return float(np.exp(np.max(self.log_atom)))

    @property
    def prior_log2_success(self) -> float:
        return float(np.max(self.log_atom) / math.log(2.0))


def log_binom_counts(n: int) -> np.ndarray:
    k = np.arange(n + 1, dtype=float)
    return gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)


def uniform_prior(bits: int) -> AggregatePrior:
    h = np.arange(bits + 1, dtype=float)
    log_count = log_binom_counts(bits)
    log_atom = np.full(bits + 1, -bits * math.log(2.0))
    log_p_h = log_count + log_atom
    return AggregatePrior(f"Uniform {bits}-bit", bits, h, log_count, log_p_h, log_atom)


def poisson_hw_prior(bits: int, lam: float) -> AggregatePrior:
    h = np.arange(bits + 1, dtype=float)
    log_count = log_binom_counts(bits)
    log_h_mass = h * math.log(lam) - gammaln(h + 1.0)
    log_p_h = log_h_mass - float(logsumexp(log_h_mass))
    log_atom = log_p_h - log_count
    return AggregatePrior(f"Poisson-HW {bits}-bit", bits, h, log_count, log_p_h, log_atom)


def gaussian_logpdf_grid(y: np.ndarray, means: np.ndarray, sigma2: float) -> np.ndarray:
    return -0.5 * math.log(2.0 * math.pi * sigma2) - ((y[:, None] - means[None, :]) ** 2) / (2.0 * sigma2)


def integration_grid(bits: int, sigma2: float, pad: float = 8.0) -> Tuple[np.ndarray, float]:
    sigma = math.sqrt(float(sigma2))
    lo = -pad * sigma
    hi = float(bits) + pad * sigma
    dx = min(0.08, max(0.01, 0.02 * sigma))
    count = int(math.ceil((hi - lo) / dx)) + 1
    y = np.linspace(lo, hi, count)
    return y, float((hi - lo) / (count - 1))


def log_integral_from_log_values(log_values: np.ndarray, dx: float) -> float:
    return float(logsumexp(log_values) + math.log(dx))


def finite_prw_log_budget(prior: AggregatePrior, sigma2: float, alpha: float) -> float:
    """README-compatible aggregate PRW budget for Gaussian noise.

    This follows the old notebook profile, using the Hamming-weight marginal
    p_H(h) inside the alpha integral.
    """
    y, dx = integration_grid(prior.bits, sigma2)
    logpdf = gaussian_logpdf_grid(y, prior.h, sigma2)
    log_inner = logsumexp(prior.log_p_h[None, :] + alpha * logpdf, axis=1)
    return alpha * log_integral_from_log_values(log_inner / alpha, dx)


def finite_prw_log2_success_bound(prior: AggregatePrior, sigma2: float, alpha: float, q_task: float) -> float:
    log_budget = finite_prw_log_budget(prior, sigma2, alpha)
    rho = invert_alpha_bound(q_task, alpha, log_budget)
    return math.log2(rho)


def normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def map_integral_uniform(bits: int, sigma2: float) -> float:
    if sigma2 <= 0.0:
        return float(bits + 1)
    sigma = math.sqrt(sigma2)
    a = 0.5 / sigma
    phi = normal_cdf(a)
    return 2.0 * phi + (bits - 1.0) * (2.0 * phi - 1.0)


def map_log2_full(prior: AggregatePrior, sigma2: float) -> float:
    if prior.label.startswith("Uniform"):
        return -float(prior.bits) + math.log2(map_integral_uniform(prior.bits, sigma2))
    y, dx = integration_grid(prior.bits, sigma2)
    logpdf = gaussian_logpdf_grid(y, prior.h, sigma2)
    log_best = np.max(prior.log_atom[None, :] + logpdf, axis=1)
    return log_integral_from_log_values(log_best, dx) / math.log(2.0)


def ball_log_counts(bits: int, radius: int) -> np.ndarray:
    out = np.full((bits + 1, bits + 1), -np.inf, dtype=float)
    log_fact = gammaln(np.arange(bits + 1, dtype=float) + 1.0)

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return -math.inf
        return float(log_fact[n] - log_fact[k] - log_fact[n - k])

    for cand_w in range(bits + 1):
        for true_w in range(bits + 1):
            vals = []
            lo = max(0, cand_w + true_w - bits)
            hi = min(cand_w, true_w)
            for overlap in range(lo, hi + 1):
                dist = cand_w + true_w - 2 * overlap
                if dist <= radius:
                    vals.append(log_choose(cand_w, overlap) + log_choose(bits - cand_w, true_w - overlap))
            if vals:
                out[cand_w, true_w] = float(logsumexp(vals))
    return out


def prior_ball_success(prior: AggregatePrior, radius: int, log_counts: np.ndarray) -> float:
    scores = logsumexp(log_counts + prior.log_atom[None, :], axis=1)
    return float(np.exp(np.max(scores)))


def map_log2_ball(prior: AggregatePrior, sigma2: float, radius: int, log_counts: np.ndarray) -> float:
    y, dx = integration_grid(prior.bits, sigma2)
    logpdf = gaussian_logpdf_grid(y, prior.h, sigma2)
    base = log_counts + prior.log_atom[None, :]
    blocks = []
    for start in range(0, len(y), 256):
        block = logpdf[start:start + 256]
        score = logsumexp(base[None, :, :] + block[:, None, :], axis=2)
        blocks.append(np.max(score, axis=1))
    return log_integral_from_log_values(np.concatenate(blocks), dx) / math.log(2.0)


def logD_alpha_bernoulli(rho: float, alpha: float, q: float) -> float:
    eps = 1e-300
    rho = min(max(float(rho), eps), 1.0 - eps)
    q = min(max(float(q), eps), 1.0 - eps)
    log_a = math.log(q) + alpha * (math.log(rho) - math.log(q))
    log_b = math.log1p(-q) + alpha * (math.log1p(-rho) - math.log1p(-q))
    return float(np.logaddexp(log_a, log_b))


def logD_alpha_bernoulli_from_logrho(log_rho: float, alpha: float, log_q: float) -> float:
    rho = math.exp(log_rho)
    q = math.exp(log_q)
    log_a = alpha * log_rho + (1.0 - alpha) * log_q
    log_b = alpha * math.log1p(-rho) + (1.0 - alpha) * math.log1p(-q)
    return float(np.logaddexp(log_a, log_b))


def invert_alpha_bound(q: float, alpha: float, log_value: float) -> float:
    if log_value <= 0.0:
        return q
    log_q = math.log(q)
    lo, hi = log_q, math.log1p(-1e-15)
    if logD_alpha_bernoulli_from_logrho(hi, alpha, log_q) <= log_value:
        return math.exp(hi)
    for _ in range(260):
        mid = 0.5 * (lo + hi)
        if logD_alpha_bernoulli_from_logrho(mid, alpha, log_q) >= log_value:
            hi = mid
        else:
            lo = mid
    return math.exp(0.5 * (lo + hi))


def required_prw_sigma2(prior: AggregatePrior, alpha: float, target_log2: float, q_task: float) -> float:
    rho = 2.0 ** target_log2
    log_target = logD_alpha_bernoulli(rho, alpha, q_task)
    return required_prw_sigma2_for_log_budget(prior, alpha, log_target)


def required_prw_sigma2_for_log_budget(prior: AggregatePrior, alpha: float, log_budget: float) -> float:
    if not math.isfinite(log_budget):
        return float("inf")
    if log_budget <= 1e-15:
        return float("inf")
    sigma2_grid, budget_grid = prw_budget_table(prior, alpha)
    # The budget decreases as sigma^2 grows.  If the small-variance endpoint is
    # already below the requested budget, zero noise is enough at this grid
    # resolution.  If the large-variance endpoint is still above it, report
    # infinity rather than extrapolating a near-prior target.
    if budget_grid[0] <= log_budget:
        return 0.0
    if budget_grid[-1] > log_budget:
        return float("inf")
    log_sigma2 = np.interp(log_budget, budget_grid[::-1], np.log(sigma2_grid)[::-1])
    return float(math.exp(log_sigma2))


def prw_budget_table(prior: AggregatePrior, alpha: float) -> Tuple[np.ndarray, np.ndarray]:
    key = (prior.label, float(alpha))
    if key in BUDGET_TABLES:
        return BUDGET_TABLES[key]
    sigma2_grid = np.geomspace(1e-4, 1e5, 240)
    budgets = np.array([finite_prw_log_budget(prior, float(s2), alpha) for s2 in sigma2_grid], dtype=float)
    # Small numerical wiggles can break interpolation; enforce monotonicity.
    for i in range(1, len(budgets)):
        if budgets[i] > budgets[i - 1]:
            budgets[i] = budgets[i - 1]
    BUDGET_TABLES[key] = (sigma2_grid, budgets)
    return BUDGET_TABLES[key]


def required_map_sigma2_uniform(bits: int, target_log2: float) -> float:
    target_integral = 2.0 ** (bits + target_log2)
    if target_integral >= bits + 1:
        return 0.0
    if target_integral <= 1.0:
        return float("inf")
    lo, hi = 0.0, 1.0
    while map_integral_uniform(bits, hi) > target_integral:
        hi *= 2.0
        if hi > 1e12:
            return float("inf")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if map_integral_uniform(bits, mid) <= target_integral:
            hi = mid
        else:
            lo = mid
    return hi


def composition_prw_sigma2(prior: AggregatePrior, alpha: float, total_target_log2: float, rounds: int, q_task: float) -> float:
    total_budget = logD_alpha_bernoulli(2.0 ** total_target_log2, alpha, q_task)
    return required_prw_sigma2_for_log_budget(prior, alpha, total_budget / rounds)


def composition_map_sigma2_uniform(bits: int, total_target_log2: float, rounds: int) -> float:
    q_log2 = -float(bits)
    local_target_log2 = q_log2 + (total_target_log2 - q_log2) / rounds
    return required_map_sigma2_uniform(bits, local_target_log2)


def font(size: int, bold: bool = False):
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
    for path in paths:
        try:
            return ImageFont.truetype(path, size)
        except Exception:
            pass
    return ImageFont.load_default()


COLORS = [
    (0, 114, 178),
    (213, 94, 0),
    (0, 158, 115),
    (204, 121, 167),
    (86, 180, 233),
    (117, 112, 179),
]


def draw_line_panel(
    d: ImageDraw.ImageDraw,
    rect: Tuple[int, int, int, int],
    title: str,
    series: Dict[str, Sequence[Tuple[float, float]]],
    x_ticks: Sequence[float],
    y_ticks: Optional[Sequence[float]],
    title_font,
    tick_font,
) -> None:
    x0, y0, x1, y1 = rect
    pw, ph = x1 - x0, y1 - y0
    all_x = [x for pts in series.values() for x, y in pts if math.isfinite(y)]
    all_y = [y for pts in series.values() for x, y in pts if math.isfinite(y)]
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if ymax <= ymin:
        ymax = ymin + 1.0
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad

    if y_ticks is None:
        y_ticks = np.linspace(ymin, ymax, 5)
    d.text((x0, y0 - 50), title, fill=(25, 25, 25), font=title_font)
    d.rectangle([x0, y0, x1, y1], fill=(250, 250, 250), outline=(172, 172, 172), width=2)

    def xm(x: float) -> float:
        return x0 + (x - xmin) / (xmax - xmin) * pw

    def ym(y: float) -> float:
        return y0 + (ymax - y) / (ymax - ymin) * ph

    for tick in x_ticks:
        px = xm(float(tick))
        d.line([(px, y0), (px, y1)], fill=(226, 226, 226), width=1)
        label = f"{tick:g}"
        tw = d.textlength(label, font=tick_font)
        d.text((px - tw / 2, y1 + 10), label, fill=(45, 45, 45), font=tick_font)

    for tick in y_ticks:
        py = ym(float(tick))
        d.line([(x0, py), (x1, py)], fill=(226, 226, 226), width=1)
        if abs(tick) >= 100:
            label = f"{tick:.0f}"
        elif abs(tick) >= 10:
            label = f"{tick:.1f}".rstrip("0").rstrip(".")
        else:
            label = f"{tick:.2g}"
        tw = d.textlength(label, font=tick_font)
        d.text((x0 - tw - 10, py - 12), label, fill=(45, 45, 45), font=tick_font)

    d.line([(x0, y1), (x1, y1)], fill=(30, 30, 30), width=3)
    d.line([(x0, y0), (x0, y1)], fill=(30, 30, 30), width=3)
    for idx, (name, pts) in enumerate(series.items()):
        color = COLORS[idx % len(COLORS)]
        clean = [(float(x), float(y)) for x, y in pts if math.isfinite(y)]
        mapped = [(xm(x), ym(y)) for x, y in clean]
        for a, b in zip(mapped, mapped[1:]):
            d.line([a, b], fill=color, width=4)
        for px, py in mapped:
            d.ellipse([px - 4, py - 4, px + 4, py + 4], fill=color, outline="white", width=1)


def draw_rotated_label(img: Image.Image, text: str, center: Tuple[int, int], text_font) -> None:
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    tw = int(measure.textlength(text, font=text_font) + 20)
    tmp = Image.new("RGBA", (tw, 64), (255, 255, 255, 0))
    td = ImageDraw.Draw(tmp)
    td.text((0, 0), text, fill=(25, 25, 25), font=text_font)
    bbox = tmp.getbbox()
    if bbox:
        tmp = tmp.crop(bbox)
    rot = tmp.rotate(90, expand=True)
    img.paste(rot, (center[0] - rot.width // 2, center[1] - rot.height // 2), rot)


def draw_legend(d: ImageDraw.ImageDraw, labels: Sequence[str], x: int, y: int, legend_font) -> None:
    cursor = x
    for idx, label in enumerate(labels):
        color = COLORS[idx % len(COLORS)]
        d.line([(cursor, y + 16), (cursor + 46, y + 16)], fill=color, width=5)
        d.text((cursor + 56, y), label, fill=(25, 25, 25), font=legend_font)
        cursor += int(66 + d.textlength(label, font=legend_font) + 34)


def log_variance_series(series: Dict[str, Sequence[Tuple[float, float]]]) -> Dict[str, List[Tuple[float, float]]]:
    transformed: Dict[str, List[Tuple[float, float]]] = {}
    for name, pts in series.items():
        transformed[name] = []
        for x, y in pts:
            if math.isfinite(float(y)):
                transformed[name].append((float(x), math.log10(1.0 + max(0.0, float(y)))))
    return transformed


def make_top_curves(rows: List[dict]) -> List[Tuple[str, Dict[str, List[Tuple[float, float]]]]]:
    sigmas = np.arange(0.25, 30.0 + 1e-9, 0.25)
    cases: List[Tuple[str, AggregatePrior, float, bool]] = []
    cases.append(("(a) Uniform 128-bit, full identification", uniform_prior(128), uniform_prior(128).prior_success, None))
    cases.append(("(b) Uniform 256-bit, full identification", uniform_prior(256), uniform_prior(256).prior_success, None))
    cases.append(("(c) Poisson-HW 256-bit, full identification", poisson_hw_prior(256, 64.0), poisson_hw_prior(256, 64.0).prior_success, None))
    prior256 = uniform_prior(256)
    ball = ball_log_counts(256, 16)
    cases.append(("(d) Uniform 256-bit, recover at least 240 bit", prior256, prior_ball_success(prior256, 16, ball), True))

    panels = []
    for title, prior, q_task, use_task_scaled_map in cases:
        curves: Dict[str, List[Tuple[float, float]]] = {f"PRW alpha={int(a)}": [] for a in CURVE_ALPHAS}
        curves["MAP / alpha=infinity"] = []
        for sigma in sigmas:
            sigma2 = float(sigma * sigma)
            for alpha in CURVE_ALPHAS:
                log2_bound = finite_prw_log2_success_bound(prior, sigma2, alpha, q_task)
                curves[f"PRW alpha={int(alpha)}"].append((float(sigma), log2_bound))
                rows.append({
                    "panel": title[:3],
                    "case": title,
                    "quantity": "posterior_log2_bound",
                    "method": f"PRW alpha={int(alpha)}",
                    "sigma": sigma,
                    "sigma2": sigma2,
                    "log2_success": log2_bound,
                    "q_task_log2": math.log2(q_task),
                })
            if use_task_scaled_map is True:
                # This mirrors the old notebook's alpha=infinity path for the
                # "miss 16" task: change q1 to the task prior success, but keep
                # the same aggregate max-density integral.
                gt = math.log2(q_task) + math.log2(map_integral_uniform(prior.bits, sigma2))
            else:
                gt = map_log2_full(prior, sigma2)
            curves["MAP / alpha=infinity"].append((float(sigma), gt))
            rows.append({
                "panel": title[:3],
                "case": title,
                "quantity": "posterior_log2_bound",
                "method": "MAP / alpha=infinity",
                "sigma": sigma,
                "sigma2": sigma2,
                "log2_success": gt,
                "q_task_log2": math.log2(q_task),
            })
        panels.append((title, curves))
    return panels


def make_bottom_curves(rows: List[dict]) -> List[Tuple[str, Dict[str, List[Tuple[float, float]]]]]:
    panels = []
    target_grid = np.arange(248.0, 252.0 + 1e-9, 0.25)
    entropy_series: Dict[str, List[Tuple[float, float]]] = {}
    for bits in [253, 254, 255, 256]:
        prior = uniform_prior(bits)
        pts = []
        for exponent in target_grid:
            sigma2 = required_prw_sigma2(prior, MAIN_ALPHA, -float(exponent), prior.prior_success)
            pts.append((float(exponent), sigma2))
            rows.append({
                "panel": "(e)",
                "case": "entropy",
                "method": f"PRW alpha={int(MAIN_ALPHA)}",
                "bits": bits,
                "target_log2": -float(exponent),
                "sigma2": sigma2,
                "quantity": "required_variance",
            })
        entropy_series[f"l={bits}"] = pts
    panels.append((r"(e) Effect from various l-bit secret entropy", entropy_series))

    prior = uniform_prior(256)
    target_grid_f = np.arange(220.0, 252.0 + 1e-9, 2.0)
    alpha_series: Dict[str, List[Tuple[float, float]]] = {}
    for alpha in CURVE_ALPHAS:
        pts = []
        for exponent in target_grid_f:
            sigma2 = required_prw_sigma2(prior, alpha, -float(exponent), prior.prior_success)
            pts.append((float(exponent), sigma2))
            rows.append({
                "panel": "(f)",
                "case": "target sweep",
                "method": f"PRW alpha={int(alpha)}",
                "bits": 256,
                "target_log2": -float(exponent),
                "sigma2": sigma2,
                "quantity": "required_variance",
            })
        alpha_series[f"PRW alpha={int(alpha)}"] = pts
    panels.append(("(f) Target sweep, 256-bit full identification", alpha_series))

    for letter, target_exponent in [("(g)", 240.0), ("(h)", 220.0)]:
        comp_series: Dict[str, List[Tuple[float, float]]] = {}
        for alpha in CURVE_ALPHAS:
            pts = []
            for t in range(1, 11):
                sigma2 = composition_prw_sigma2(prior, alpha, -target_exponent, t, prior.prior_success)
                pts.append((float(t), sigma2))
                rows.append({
                    "panel": letter,
                    "case": "composition",
                    "method": f"PRW alpha={int(alpha)}",
                    "bits": 256,
                    "total_target_log2": -target_exponent,
                    "rounds": t,
                    "sigma2": sigma2,
                    "quantity": "per_round_required_variance",
                })
            comp_series[f"PRW alpha={int(alpha)}"] = pts
        panels.append((fr"{letter} Composition, full identification 2^-{int(target_exponent)}", comp_series))
    return panels


def draw_combined_2x4(path: Path, top_panels, bottom_panels) -> None:
    width, height = 5600, 2600
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    title_font = font(46, True)
    tick_font = font(34)
    axis_font = font(48, True)
    legend_font = font(42)
    margin_left, margin_right = 210, 70
    gap_x, gap_y = 125, 260
    panel_w = (width - margin_left - margin_right - 3 * gap_x) // 4
    panel_h = 820
    top_y, bottom_y = 270, 1580

    draw_legend(d, list(top_panels[0][1].keys()), 1150, 60, legend_font)
    draw_legend(d, list(bottom_panels[2][1].keys()), 1650, 1375, legend_font)

    for row_idx, panels in enumerate([top_panels, bottom_panels]):
        y = top_y if row_idx == 0 else bottom_y
        for col, (title, series) in enumerate(panels):
            x = margin_left + col * (panel_w + gap_x)
            plot_series = series if row_idx == 0 else log_variance_series(series)
            if row_idx == 0:
                xticks = [0, 10, 20, 30]
                yticks = None
            elif col < 2:
                xs = [x for pts in plot_series.values() for x, _ in pts]
                xticks = np.linspace(min(xs), max(xs), 5)
                yticks = None
            else:
                xticks = list(range(1, 11))
                yticks = None
            draw_line_panel(d, (x, y, x + panel_w, y + panel_h), title, plot_series, xticks, yticks, title_font, tick_font)
            if row_idx == 1 and col == 0:
                label_x = x + 26
                label_y = y + 26
                for idx, name in enumerate(series.keys()):
                    color = COLORS[idx % len(COLORS)]
                    yy = label_y + idx * 42
                    d.line([(label_x, yy + 15), (label_x + 44, yy + 15)], fill=color, width=5)
                    d.text((label_x + 54, yy), name, fill=(25, 25, 25), font=legend_font)

    d.text((width // 2 - int(d.textlength("Gaussian noise std sigma", font=axis_font)) // 2, top_y + panel_h + 90), "Gaussian noise std sigma", fill=(25, 25, 25), font=axis_font)
    d.text((margin_left + panel_w + gap_x // 2, bottom_y + panel_h + 90), "log2 of target posterior success rate", fill=(25, 25, 25), font=axis_font)
    d.text((margin_left + 2 * (panel_w + gap_x) + panel_w // 2, bottom_y + panel_h + 90), "Composition iteration", fill=(25, 25, 25), font=axis_font)
    draw_rotated_label(img, "log2 of posterior success rate", (80, top_y + panel_h // 2), axis_font)
    draw_rotated_label(img, "log10(1 + Gaussian noise variance sigma^2)", (80, bottom_y + panel_h // 2), axis_font)

    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def draw_variance_comparison(rows: List[dict]) -> None:
    prior = uniform_prior(256)
    target_exponents = np.arange(220.0, 252.0 + 1e-9, 2.0)
    left_series: Dict[str, List[Tuple[float, float]]] = {}
    for alpha in CURVE_ALPHAS:
        pts = []
        for exponent in target_exponents:
            sigma2 = required_prw_sigma2(prior, alpha, -float(exponent), prior.prior_success)
            pts.append((float(exponent), sigma2))
            rows.append({
                "panel": "variance comparison one-shot",
                "method": f"PRW alpha={int(alpha)}",
                "target_log2": -float(exponent),
                "sigma2": sigma2,
            })
        left_series[f"PRW alpha={int(alpha)}"] = pts
    map_pts = []
    for exponent in target_exponents:
        sigma2 = required_map_sigma2_uniform(256, -float(exponent))
        map_pts.append((float(exponent), sigma2))
        rows.append({
            "panel": "variance comparison one-shot",
            "method": "MAP / alpha=infinity",
            "target_log2": -float(exponent),
            "sigma2": sigma2,
        })
    left_series["MAP / alpha=infinity"] = map_pts

    right_series: Dict[str, List[Tuple[float, float]]] = {}
    for alpha in CURVE_ALPHAS:
        pts = []
        for t in range(1, 11):
            sigma2 = composition_prw_sigma2(prior, alpha, -240.0, t, prior.prior_success)
            pts.append((float(t), sigma2))
            rows.append({
                "panel": "variance comparison composition",
                "method": f"PRW alpha={int(alpha)}",
                "total_target_log2": -240.0,
                "rounds": t,
                "sigma2": sigma2,
            })
        right_series[f"PRW alpha={int(alpha)}"] = pts
    map_comp = []
    for t in range(1, 11):
        sigma2 = composition_map_sigma2_uniform(256, -240.0, t)
        map_comp.append((float(t), sigma2))
        rows.append({
            "panel": "variance comparison composition",
            "method": "MAP / alpha=infinity",
            "total_target_log2": -240.0,
            "rounds": t,
            "sigma2": sigma2,
        })
    right_series["MAP / alpha=infinity"] = map_comp

    width, height = 2800, 1150
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    title_font = font(46, True)
    tick_font = font(34)
    axis_font = font(46, True)
    legend_font = font(40)
    draw_legend(d, list(left_series.keys()), 420, 55, legend_font)
    draw_line_panel(
        d,
        (210, 250, 1320, 930),
        "(i) One-shot variance, 256-bit full identification",
        log_variance_series(left_series),
        [220, 228, 236, 244, 252],
        None,
        title_font,
        tick_font,
    )
    draw_line_panel(
        d,
        (1570, 250, 2680, 930),
        "(j) Composition variance, target 2^-240",
        log_variance_series(right_series),
        list(range(1, 11)),
        None,
        title_font,
        tick_font,
    )
    d.text((530, 1010), "log2 of target posterior success rate", fill=(25, 25, 25), font=axis_font)
    d.text((1850, 1010), "Composition iteration", fill=(25, 25, 25), font=axis_font)
    draw_rotated_label(img, "log10(1 + Gaussian noise variance sigma^2)", (80, 585), axis_font)
    img.save(OUT / "prw_vs_map_variance.png")


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


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    top = make_top_curves(rows)
    bottom = make_bottom_curves(rows)
    draw_combined_2x4(OUT / "power_prw_finite_alpha_2x4.png", top, bottom)
    draw_variance_comparison(rows)
    write_csv(OUT / "power_prw_finite_alpha_results.csv", rows)
    print(f"wrote {OUT / 'power_prw_finite_alpha_2x4.png'}")
    print(f"wrote {OUT / 'prw_vs_map_variance.png'}")
    print(f"wrote {OUT / 'power_prw_finite_alpha_results.csv'}")


if __name__ == "__main__":
    run()
