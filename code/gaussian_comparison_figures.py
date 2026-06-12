"""Gaussian-noise comparison figures for aggregate Hamming-weight leakage.

The first four figures compare:
  1. PAC-KL/MI,
  2. PAC alpha divergence with the best alpha selected per noise level,
  3. our PRW bound at fixed alpha values,
  4. exact ground truth.

All curves use independent continuous Gaussian noise added to the scalar
aggregate Hamming-weight leakage.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.special import gammaln, logsumexp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "gaussian_comparison"
PAC_INVERSION_STEPS = 240


@dataclass
class AggregatePrior:
    name: str
    bits: int
    weights: np.ndarray
    log_count: np.ndarray
    log_atom: np.ndarray
    log_p_h: np.ndarray

    @property
    def prior_identification_success(self) -> float:
        return float(np.exp(np.max(self.log_atom)))

    @property
    def mean(self) -> float:
        return float(np.sum(self.weights * np.exp(self.log_p_h)))


def log_binom_counts(n: int) -> np.ndarray:
    k = np.arange(n + 1, dtype=float)
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def uniform_prior(bits: int) -> AggregatePrior:
    weights = np.arange(bits + 1, dtype=float)
    lc = log_binom_counts(bits)
    log_atom = np.full(bits + 1, -bits * math.log(2.0), dtype=float)
    log_p_h = lc + log_atom
    return AggregatePrior(f"uniform {bits}-bit", bits, weights, lc, log_atom, log_p_h)


def hamming_discrete_gaussian_prior(bits: int, tau: float, center: Optional[float] = None) -> AggregatePrior:
    """Per-string prior tilted by a discrete Gaussian in Hamming weight."""
    weights = np.arange(bits + 1, dtype=float)
    lc = log_binom_counts(bits)
    c = bits / 2.0 if center is None else float(center)
    tilt = -((weights - c) ** 2) / (2.0 * tau * tau)
    log_z = float(logsumexp(lc + tilt))
    log_atom = tilt - log_z
    log_p_h = lc + log_atom
    return AggregatePrior(f"discrete Gaussian center={c:g}, tau={tau:g}", bits, weights, lc, log_atom, log_p_h)


def hamming_poisson_prior(bits: int, lam: float) -> AggregatePrior:
    """Hamming weight H follows a truncated Poisson law; x is uniform given H."""
    if lam <= 0.0:
        raise ValueError("Poisson mean must be positive")
    weights = np.arange(bits + 1, dtype=float)
    lc = log_binom_counts(bits)
    log_h_mass = weights * math.log(lam) - gammaln(weights + 1.0)
    log_p_h = log_h_mass - float(logsumexp(log_h_mass))
    log_atom = log_p_h - lc
    return AggregatePrior(f"Poisson-HW lambda={lam:g}", bits, weights, lc, log_atom, log_p_h)


def gaussian_logpdf_grid(y: np.ndarray, means: np.ndarray, sigma: float) -> np.ndarray:
    return -0.5 * math.log(2.0 * math.pi * sigma * sigma) - ((y[:, None] - means[None, :]) ** 2) / (2.0 * sigma * sigma)


def integration_grid(bits: int, sigma: float, pad: float = 8.0) -> Tuple[np.ndarray, float]:
    lo = -pad * sigma
    hi = bits + pad * sigma
    dx = min(0.20, max(0.05, 0.04 * sigma))
    count = int(math.ceil((hi - lo) / dx)) + 1
    return np.linspace(lo, hi, count), (hi - lo) / (count - 1)


def log_integral_from_log_values(log_values: np.ndarray, dx: float) -> float:
    return float(logsumexp(log_values) + math.log(dx))


def log_ground_truth_full(prior: AggregatePrior, sigma: float) -> float:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    log_best = np.max(prior.log_atom[None, :] + logpdf, axis=1)
    return log_integral_from_log_values(log_best, dx)


def ball_log_counts(bits: int, radius: int) -> np.ndarray:
    log_counts = np.full((bits + 1, bits + 1), -np.inf, dtype=float)
    log_fact = gammaln(np.arange(bits + 1, dtype=float) + 1.0)

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return -math.inf
        return float(log_fact[n] - log_fact[k] - log_fact[n - k])

    for cand_weight in range(bits + 1):
        for true_weight in range(bits + 1):
            vals = []
            lo = max(0, cand_weight + true_weight - bits)
            hi = min(cand_weight, true_weight)
            for overlap in range(lo, hi + 1):
                dist = cand_weight + true_weight - 2 * overlap
                if dist <= radius:
                    vals.append(log_choose(cand_weight, overlap) + log_choose(bits - cand_weight, true_weight - overlap))
            if vals:
                log_counts[cand_weight, true_weight] = float(logsumexp(vals))
    return log_counts


def log_prior_task_success(prior: AggregatePrior, radius: int) -> float:
    log_counts = ball_log_counts(prior.bits, radius)
    scores = logsumexp(log_counts + prior.log_atom[None, :], axis=1)
    return float(np.max(scores))


def log_ground_truth_ball(prior: AggregatePrior, sigma: float, radius: int, precomputed_log_counts: np.ndarray) -> float:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    # score_j(y) = sum_k count(j,k) atom(k) phi(y-k)
    log_score = []
    base = precomputed_log_counts + prior.log_atom[None, :]
    for start in range(0, len(y), 512):
        block = logpdf[start:start + 512]
        vals = logsumexp(base[None, :, :] + block[:, None, :], axis=2)
        log_score.append(np.max(vals, axis=1))
    return log_integral_from_log_values(np.concatenate(log_score), dx)


def ball_ground_truth_and_prw(
    prior: AggregatePrior,
    sigma: float,
    alpha_values: Sequence[float],
    precomputed_log_counts: np.ndarray,
) -> Tuple[float, Dict[float, float]]:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    base = precomputed_log_counts + prior.log_atom[None, :]
    gt_blocks = []
    prw_blocks: Dict[float, List[np.ndarray]] = {float(alpha): [] for alpha in alpha_values}
    for start in range(0, len(y), 128):
        block = logpdf[start:start + 128]
        log_score = logsumexp(base[None, :, :] + block[:, None, :], axis=2)
        gt_blocks.append(np.max(log_score, axis=1))
        for alpha in prw_blocks:
            log_norm = logsumexp(prior.log_count[None, :] + alpha * log_score, axis=1) / alpha
            prw_blocks[alpha].append(log_norm)
    gt = log_integral_from_log_values(np.concatenate(gt_blocks), dx)
    prw = {
        alpha: log_integral_from_log_values(np.concatenate(blocks), dx)
        for alpha, blocks in prw_blocks.items()
    }
    return gt, prw


def ball_prw_only(
    prior: AggregatePrior,
    sigma: float,
    alpha_values: Sequence[float],
    precomputed_log_counts: np.ndarray,
) -> Dict[float, float]:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    base = precomputed_log_counts + prior.log_atom[None, :]
    prw_blocks: Dict[float, List[np.ndarray]] = {float(alpha): [] for alpha in alpha_values}
    for start in range(0, len(y), 128):
        block = logpdf[start:start + 128]
        log_score = logsumexp(base[None, :, :] + block[:, None, :], axis=2)
        for alpha in prw_blocks:
            log_norm = logsumexp(prior.log_count[None, :] + alpha * log_score, axis=1) / alpha
            prw_blocks[alpha].append(log_norm)
    return {
        alpha: log_integral_from_log_values(np.concatenate(blocks), dx)
        for alpha, blocks in prw_blocks.items()
    }


def mutual_information(prior: AggregatePrior, sigma: float) -> float:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    log_joint = prior.log_p_h[None, :] + logpdf
    log_py = logsumexp(log_joint, axis=1)
    integrand = np.sum(np.exp(log_joint) * (logpdf - log_py[:, None]), axis=1)
    return float(np.sum(integrand) * dx)


def log_pac_alpha_cost(prior: AggregatePrior, sigma: float, alpha: float) -> float:
    mu = float(np.sum(prior.weights * np.exp(prior.log_p_h)))
    c = alpha * (alpha - 1.0) / 2.0
    return float(logsumexp(prior.log_p_h + c * ((prior.weights - mu) ** 2) / (sigma * sigma)))


def log_pac_alpha_cost_vector(prior: AggregatePrior, sigma: float, alpha_values: np.ndarray) -> np.ndarray:
    mu = float(np.sum(prior.weights * np.exp(prior.log_p_h)))
    scaled_sq = ((prior.weights - mu) ** 2) / (sigma * sigma)
    c = alpha_values * (alpha_values - 1.0) / 2.0
    return logsumexp(prior.log_p_h[None, :] + c[:, None] * scaled_sq[None, :], axis=1)


def log_prw_alpha_cost(prior: AggregatePrior, sigma: float, alpha: float) -> float:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    log_inner = logsumexp(prior.log_count[None, :] + alpha * prior.log_atom[None, :] + alpha * logpdf, axis=1)
    log_i = log_integral_from_log_values(log_inner / alpha, dx)
    return alpha * log_i


def log_prw_alpha_ball(
    prior: AggregatePrior,
    sigma: float,
    radius: int,
    alpha: float,
    precomputed_log_counts: np.ndarray,
) -> float:
    y, dx = integration_grid(prior.bits, sigma)
    logpdf = gaussian_logpdf_grid(y, prior.weights, sigma)
    base = precomputed_log_counts + prior.log_atom[None, :]
    log_norm_blocks = []
    for start in range(0, len(y), 384):
        block = logpdf[start:start + 384]
        log_score = logsumexp(base[None, :, :] + block[:, None, :], axis=2)
        log_norm = logsumexp(prior.log_count[None, :] + alpha * log_score, axis=1) / alpha
        log_norm_blocks.append(log_norm)
    return log_integral_from_log_values(np.concatenate(log_norm_blocks), dx)


def logD_alpha_bernoulli(rho: float, alpha: float, q: float) -> float:
    eps = 1e-300
    rho = min(max(float(rho), eps), 1.0 - eps)
    q = min(max(float(q), eps), 1.0 - eps)
    log_a = math.log(q) + alpha * (math.log(rho) - math.log(q))
    log_b = math.log1p(-q) + alpha * (math.log1p(-rho) - math.log1p(-q))
    return float(logsumexp([log_a, log_b]))


def logD_alpha_bernoulli_from_logrho(log_rho: float, alpha: float, log_q: float) -> float:
    rho = math.exp(log_rho)
    q = math.exp(log_q)
    log_a = alpha * log_rho + (1.0 - alpha) * log_q
    log_b = alpha * math.log1p(-rho) + (1.0 - alpha) * math.log1p(-q)
    return float(logsumexp([log_a, log_b]))


def binary_kl(rho: float, q: float) -> float:
    eps = 1e-300
    rho = min(max(float(rho), eps), 1.0 - eps)
    q = min(max(float(q), eps), 1.0 - eps)
    return rho * math.log(rho / q) + (1.0 - rho) * math.log((1.0 - rho) / (1.0 - q))


def binary_kl_from_logrho(log_rho: float, log_q: float) -> float:
    rho = math.exp(log_rho)
    q = math.exp(log_q)
    return rho * (log_rho - log_q) + (1.0 - rho) * (math.log1p(-rho) - math.log1p(-q))


def invert_kl_bound(q: float, value: float) -> float:
    if value <= 0.0:
        return q
    log_q = math.log(q)
    lo, hi = log_q, math.log1p(-1e-15)
    if binary_kl_from_logrho(hi, log_q) <= value:
        return math.exp(hi)
    for _ in range(PAC_INVERSION_STEPS):
        mid = 0.5 * (lo + hi)
        if binary_kl_from_logrho(mid, log_q) >= value:
            hi = mid
        else:
            lo = mid
    return math.exp(0.5 * (lo + hi))


def invert_alpha_bound(q: float, alpha: float, log_value: float) -> float:
    if log_value <= 0.0:
        return q
    log_q = math.log(q)
    lo, hi = log_q, math.log1p(-1e-15)
    if logD_alpha_bernoulli_from_logrho(hi, alpha, log_q) <= log_value:
        return math.exp(hi)
    for _ in range(PAC_INVERSION_STEPS):
        mid = 0.5 * (lo + hi)
        if logD_alpha_bernoulli_from_logrho(mid, alpha, log_q) >= log_value:
            hi = mid
        else:
            lo = mid
    return math.exp(0.5 * (lo + hi))


def invert_alpha_bound_logrho_vector(q: float, alpha_values: np.ndarray, log_values: np.ndarray) -> np.ndarray:
    log_q = math.log(q)
    lo = np.full(alpha_values.shape, log_q, dtype=float)
    hi = np.full(alpha_values.shape, math.log1p(-1e-15), dtype=float)
    log_one_minus_q = math.log1p(-q)

    for _ in range(PAC_INVERSION_STEPS):
        mid = 0.5 * (lo + hi)
        rho = np.exp(mid)
        log_a = alpha_values * mid + (1.0 - alpha_values) * log_q
        log_b = alpha_values * np.log1p(-rho) + (1.0 - alpha_values) * log_one_minus_q
        log_d = np.logaddexp(log_a, log_b)
        mask = log_d >= log_values
        hi = np.where(mask, mid, hi)
        lo = np.where(mask, lo, mid)
    return 0.5 * (lo + hi)


def pac_alpha_best_log2(prior: AggregatePrior, sigma: float, q_task: float, alpha_grid: Sequence[float]) -> Tuple[float, float]:
    best_log2 = 0.0
    best_alpha = float(np.asarray(alpha_grid, dtype=float)[0])

    def evaluate_many(alpha_values: np.ndarray) -> None:
        nonlocal best_log2, best_alpha
        alpha_values = np.unique(np.asarray(alpha_values, dtype=float))
        alpha_values = alpha_values[alpha_values > 1.0]
        if alpha_values.size == 0:
            return
        log_values = log_pac_alpha_cost_vector(prior, sigma, alpha_values)
        log_rhos = invert_alpha_bound_logrho_vector(q_task, alpha_values, log_values)
        vals = log_rhos / math.log(2.0)
        idx = int(np.argmin(vals))
        val = float(vals[idx])
        if val < best_log2:
            best_log2 = val
            best_alpha = float(alpha_values[idx])

    evaluate_many(np.asarray(alpha_grid, dtype=float))
    for width, count in ((0.40, 101), (0.08, 81), (0.02, 61)):
        lo = max(1.001, best_alpha - width)
        hi = best_alpha + width
        evaluate_many(np.linspace(lo, hi, count))
    return best_log2, best_alpha


def font(size: int, bold: bool = False):
    paths = []
    if bold:
        paths.extend(["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/Library/Fonts/Arial Bold.ttf"])
    paths.extend(["/System/Library/Fonts/Supplemental/Arial.ttf", "/Library/Fonts/Arial.ttf"])
    for p in paths:
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            pass
    return ImageFont.load_default()


def plot_curves(path: Path, title: str, subtitle: str, sigmas: np.ndarray, curves: Dict[str, Sequence[float]]) -> None:
    width, height = 1500, 930
    margin = {"left": 150, "right": 360, "top": 120, "bottom": 125}
    pw, ph = width - margin["left"] - margin["right"], height - margin["top"] - margin["bottom"]
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    ft, fs, fa, fk = font(32, True), font(18), font(23), font(18)
    d.text((margin["left"], 34), title, fill=(25, 25, 25), font=ft)
    d.text((margin["left"], 78), subtitle, fill=(80, 80, 80), font=fs)

    xmin, xmax = float(np.min(sigmas)), float(np.max(sigmas))
    all_y = [float(v) for vals in curves.values() for v in vals if math.isfinite(float(v))]
    ymin = math.floor(min(all_y) - 1.0)
    ymax = math.ceil(max(all_y) + 1.0)

    def xm(x: float) -> float:
        return margin["left"] + (x - xmin) / (xmax - xmin) * pw

    def ym(y: float) -> float:
        return margin["top"] + (ymax - y) / (ymax - ymin) * ph

    d.rectangle([margin["left"], margin["top"], margin["left"] + pw, margin["top"] + ph], fill=(250, 250, 250), outline=(180, 180, 180))
    for x in np.linspace(xmin, xmax, 6):
        px = xm(float(x))
        d.line([(px, margin["top"]), (px, margin["top"] + ph)], fill=(226, 226, 226))
        label = f"{x:.0f}"
        tw = d.textlength(label, font=fk)
        d.text((px - tw / 2, margin["top"] + ph + 14), label, fill=(45, 45, 45), font=fk)
    y_ticks = np.linspace(ymin, ymax, 7)
    for y in y_ticks:
        py = ym(float(y))
        d.line([(margin["left"], py), (margin["left"] + pw, py)], fill=(226, 226, 226))
        label = f"{y:.0f}"
        tw = d.textlength(label, font=fk)
        d.text((margin["left"] - tw - 12, py - 9), label, fill=(45, 45, 45), font=fk)
    d.line([(margin["left"], margin["top"] + ph), (margin["left"] + pw, margin["top"] + ph)], fill=(30, 30, 30), width=2)
    d.line([(margin["left"], margin["top"]), (margin["left"], margin["top"] + ph)], fill=(30, 30, 30), width=2)

    colors = [
        (0, 114, 178),
        (213, 94, 0),
        (0, 158, 115),
        (230, 159, 0),
        (204, 121, 167),
        (150, 150, 150),
    ]
    curve_items = list(curves.items())
    draw_order = list(range(len(curve_items)))
    for i, (name, _) in enumerate(curve_items):
        if name == "Ground truth":
            draw_order = [i] + [j for j in draw_order if j != i]
            break
    for idx in draw_order:
        name, vals = curve_items[idx]
        color = colors[idx % len(colors)]
        pts = [(xm(float(x)), ym(float(y))) for x, y in zip(sigmas, vals) if math.isfinite(float(y))]
        for a, b in zip(pts, pts[1:]):
            d.line([a, b], fill=color, width=4)
        for px, py in pts[:: max(1, len(pts) // 12)]:
            d.ellipse([px - 4, py - 4, px + 4, py + 4], fill=color, outline="white", width=1)
        lx = margin["left"] + pw + 35
        ly = margin["top"] + 22 + idx * 36
        d.line([(lx, ly + 10), (lx + 48, ly + 10)], fill=color, width=4)
        d.text((lx + 60, ly), name, fill=(35, 35, 35), font=fk)

    xlabel = "Gaussian noise std σ"
    ylabel = "Logarithm log_2 posterior success"
    tw = d.textlength(xlabel, font=fa)
    d.text((margin["left"] + pw / 2 - tw / 2, height - 74), xlabel, fill=(30, 30, 30), font=fa)
    label_img = Image.new("RGBA", (760, 40), (255, 255, 255, 0))
    ld = ImageDraw.Draw(label_img)
    ld.text((0, 0), ylabel, fill=(30, 30, 30), font=fa)
    rot = label_img.rotate(90, expand=True)
    img.paste(rot, (30, margin["top"] + ph // 2 - rot.height // 2), rot)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def display_curve_label(name: str) -> str:
    return name.replace("alpha", "α")


def plot_combined_curves(
    path: Path,
    panels: Sequence[tuple],
) -> None:
    final_width, final_height = 5600, 980
    render_scale = 2
    width, height = final_width * render_scale, final_height * render_scale
    margin = {
        "left": 185 * render_scale,
        "right": 60 * render_scale,
        "top": 170 * render_scale,
        "bottom": 125 * render_scale,
    }
    gap_x = 95 * render_scale
    panel_w = (width - margin["left"] - margin["right"] - 3 * gap_x) / 4.0
    panel_h = height - margin["top"] - margin["bottom"]

    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)
    title_font = font(40 * render_scale, True)
    tick_font = font(33 * render_scale)
    x_axis_font = font(50 * render_scale, True)
    y_axis_font = x_axis_font
    legend_font = font(38 * render_scale)

    colors = [
        (0, 114, 178),
        (213, 94, 0),
        (0, 158, 115),
        (230, 159, 0),
        (204, 121, 167),
        (150, 150, 150),
    ]
    curve_names = list(panels[0][2].keys())

    def draw_zoom_inset(
        panel_box: Tuple[float, float, float, float],
        zoom_sigmas: np.ndarray,
        zoom_curves: Dict[str, Sequence[float]],
        avoid_polylines: Sequence[Sequence[Tuple[float, float]]],
    ) -> None:
        inset_names = [
            "PAC-alpha(best)",
            "PRW alpha=16",
            "PRW alpha=64",
            "PRW alpha=256",
        ]
        inset_names = [name for name in inset_names if name in zoom_curves]
        if not inset_names:
            return

        px0, py0, px1, py1 = panel_box
        pw = px1 - px0
        ph = py1 - py0

        def orient(a: Tuple[float, float], b: Tuple[float, float], c: Tuple[float, float]) -> float:
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

        def segments_intersect(
            a: Tuple[float, float],
            b: Tuple[float, float],
            c: Tuple[float, float],
            e: Tuple[float, float],
        ) -> bool:
            o1 = orient(a, b, c)
            o2 = orient(a, b, e)
            o3 = orient(c, e, a)
            o4 = orient(c, e, b)
            return (o1 == 0.0 or o2 == 0.0 or o1 * o2 < 0.0) and (o3 == 0.0 or o4 == 0.0 or o3 * o4 < 0.0)

        def segment_hits_rect(a: Tuple[float, float], b: Tuple[float, float], rect: Tuple[float, float, float, float]) -> bool:
            left, top, right, bottom = rect
            if max(a[0], b[0]) < left or min(a[0], b[0]) > right or max(a[1], b[1]) < top or min(a[1], b[1]) > bottom:
                return False
            if left <= a[0] <= right and top <= a[1] <= bottom:
                return True
            if left <= b[0] <= right and top <= b[1] <= bottom:
                return True
            corners = [(left, top), (right, top), (right, bottom), (left, bottom)]
            edges = list(zip(corners, corners[1:] + corners[:1]))
            return any(segments_intersect(a, b, edge_a, edge_b) for edge_a, edge_b in edges)

        def rect_hits_curves(rect: Tuple[float, float, float, float], pad: float) -> bool:
            padded = (rect[0] - pad, rect[1] - pad, rect[2] + pad, rect[3] + pad)
            for pts in avoid_polylines:
                for a, b in zip(pts, pts[1:]):
                    if segment_hits_rect(a, b, padded):
                        return True
            return False

        def choose_rect() -> Tuple[float, float, float, float]:
            preferred = (px0 + 0.70 * pw, py0 + 0.46 * ph)
            best: Optional[Tuple[float, float, float, float]] = None
            best_score = float("inf")
            for factor in [1.0, 0.94, 0.88, 0.82, 0.76, 0.70]:
                iw = 0.60 * pw * factor
                ih = 0.34 * ph * factor
                x_candidates = np.linspace(px0 + 0.04 * pw, px1 - iw - 0.02 * pw, 10)
                y_candidates = np.linspace(py0 + 0.08 * ph, py1 - ih - 0.05 * ph, 14)
                for ix0_candidate in x_candidates:
                    for iy0_candidate in y_candidates:
                        rect = (float(ix0_candidate), float(iy0_candidate), float(ix0_candidate + iw), float(iy0_candidate + ih))
                        if rect_hits_curves(rect, 9 * render_scale):
                            continue
                        center = ((rect[0] + rect[2]) / 2.0, (rect[1] + rect[3]) / 2.0)
                        score = abs(center[0] - preferred[0]) + 0.75 * abs(center[1] - preferred[1])
                        if score < best_score:
                            best_score = score
                            best = rect
                if best is not None:
                    return best
            return (px1 - 0.42 * pw, py0 + 0.30 * ph, px1 - 0.02 * pw, py0 + 0.56 * ph)

        ix0, iy0, ix1, iy1 = choose_rect()
        iw = ix1 - ix0
        ih = iy1 - iy0

        zxmin = float(np.min(zoom_sigmas))
        zxmax = float(np.max(zoom_sigmas))
        y_values = [
            float(value)
            for name in inset_names
            for value in zoom_curves[name]
            if math.isfinite(float(value))
        ]
        if not y_values:
            return
        zymin = min(y_values)
        zymax = max(y_values)
        y_pad = max(0.4, 0.08 * (zymax - zymin))
        zymin = math.floor(zymin - y_pad)
        zymax = math.ceil(zymax + y_pad)
        if zymax <= zymin:
            zymax = zymin + 1.0

        def zxm(x: float) -> float:
            return ix0 + (x - zxmin) / (zxmax - zxmin) * iw

        def zym(y: float) -> float:
            return iy0 + (zymax - y) / (zymax - zymin) * ih

        inset_tick_font = font(24 * render_scale)
        inset_title_font = font(25 * render_scale, True)
        d.rectangle([ix0, iy0, ix1, iy1], fill=(255, 255, 255), outline=(120, 120, 120), width=2 * render_scale)
        for x in [40, 60, 80]:
            px = zxm(float(x))
            d.line([(px, iy0), (px, iy1)], fill=(232, 232, 232), width=1 * render_scale)
            label = f"{x:g}"
            tw = d.textlength(label, font=inset_tick_font)
            d.text((px - tw / 2, iy1 - 31 * render_scale), label, fill=(55, 55, 55), font=inset_tick_font)
        for y in np.linspace(zymin, zymax, 3):
            py = zym(float(y))
            d.line([(ix0, py), (ix1, py)], fill=(232, 232, 232), width=1 * render_scale)
            label = f"{y:.0f}"
            d.text((ix0 + 8 * render_scale, py - 13 * render_scale), label, fill=(55, 55, 55), font=inset_tick_font)

        d.text((ix0 + 10 * render_scale, iy0 + 7 * render_scale), "σ=40-80", fill=(35, 35, 35), font=inset_title_font)
        for name in inset_names:
            idx = curve_names.index(name)
            color = colors[idx % len(colors)]
            pts = [
                (zxm(float(x)), zym(float(y)))
                for x, y in zip(zoom_sigmas, zoom_curves[name])
                if math.isfinite(float(y))
            ]
            for a, b in zip(pts, pts[1:]):
                d.line([a, b], fill=color, width=3 * render_scale)

    legend_items = []
    for idx, name in enumerate(curve_names):
        label = display_curve_label(name)
        legend_items.append((label, colors[idx % len(colors)]))
    legend_line = 50 * render_scale
    legend_gap = 34 * render_scale
    legend_text_gap = 14 * render_scale
    total_legend_width = sum(legend_line + legend_text_gap + d.textlength(label, font=legend_font) for label, _ in legend_items) + legend_gap * (len(legend_items) - 1)
    lx = (width - total_legend_width) / 2.0
    ly = 38 * render_scale
    for label, color in legend_items:
        d.line([(lx, ly + 20 * render_scale), (lx + legend_line, ly + 20 * render_scale)], fill=color, width=7 * render_scale)
        d.text((lx + legend_line + legend_text_gap, ly), label, fill=(30, 30, 30), font=legend_font)
        lx += legend_line + legend_text_gap + d.textlength(label, font=legend_font) + legend_gap

    x_ticks = [0.25, 10, 20, 30]
    for panel_idx, panel in enumerate(panels):
        title, sigmas, curves = panel[:3]
        zoom_sigmas = panel[3] if len(panel) >= 5 else None
        zoom_curves = panel[4] if len(panel) >= 5 else None
        col = panel_idx
        x0 = margin["left"] + col * (panel_w + gap_x)
        y0 = margin["top"]
        x1 = x0 + panel_w
        y1 = y0 + panel_h
        xmin, xmax = float(np.min(sigmas)), float(np.max(sigmas))
        all_y = [float(v) for vals in curves.values() for v in vals if math.isfinite(float(v))]
        ymin = math.floor(min(all_y) - 1.0)
        ymax = math.ceil(max(all_y) + 1.0)
        if ymax <= ymin:
            ymax = ymin + 1

        def xm(x: float) -> float:
            return x0 + (x - xmin) / (xmax - xmin) * panel_w

        def ym(y: float) -> float:
            return y0 + (ymax - y) / (ymax - ymin) * panel_h

        d.text((x0, y0 - 58 * render_scale), title, fill=(25, 25, 25), font=title_font)
        d.rectangle([x0, y0, x1, y1], fill=(250, 250, 250), outline=(178, 178, 178))

        for x in x_ticks:
            px = xm(float(x))
            d.line([(px, y0), (px, y1)], fill=(226, 226, 226))
            label = f"{x:g}"
            tw = d.textlength(label, font=tick_font)
            d.text((px - tw / 2, y1 + 12 * render_scale), label, fill=(45, 45, 45), font=tick_font)

        y_ticks = np.linspace(ymin, ymax, 6)
        for y in y_ticks:
            py = ym(float(y))
            d.line([(x0, py), (x1, py)], fill=(226, 226, 226))
            label = f"{y:.0f}"
            tw = d.textlength(label, font=tick_font)
            d.text((x0 - tw - 13 * render_scale, py - 13 * render_scale), label, fill=(45, 45, 45), font=tick_font)

        d.line([(x0, y1), (x1, y1)], fill=(30, 30, 30), width=3 * render_scale)
        d.line([(x0, y0), (x0, y1)], fill=(30, 30, 30), width=3 * render_scale)

        curve_items = list(curves.items())
        draw_order = list(range(len(curve_items)))
        for i, (name, _) in enumerate(curve_items):
            if name == "Ground truth":
                draw_order = [i] + [j for j in draw_order if j != i]
                break
        for idx in draw_order:
            name, vals = curve_items[idx]
            color = colors[idx % len(colors)]
            pts = [(xm(float(x)), ym(float(y))) for x, y in zip(sigmas, vals) if math.isfinite(float(y))]
            for a, b in zip(pts, pts[1:]):
                d.line([a, b], fill=color, width=5 * render_scale)

        if zoom_sigmas is not None and zoom_curves is not None:
            avoid_polylines = [
                [(xm(float(x)), ym(float(y))) for x, y in zip(sigmas, vals) if math.isfinite(float(y))]
                for _, vals in curve_items
            ]
            draw_zoom_inset((x0, y0, x1, y1), zoom_sigmas, zoom_curves, avoid_polylines)

    xlabel = "Gaussian noise std σ"
    tw = d.textlength(xlabel, font=x_axis_font)
    d.text((width / 2 - tw / 2, height - 66 * render_scale), xlabel, fill=(25, 25, 25), font=x_axis_font)

    prefix = "log"
    suffix = " of posterior success rate"
    sub_font = font(39 * render_scale, True)
    prefix_w = d.textlength(prefix, font=y_axis_font)
    sub_w = d.textlength("2", font=sub_font)
    suffix_w = d.textlength(suffix, font=y_axis_font)
    label_img = Image.new("RGBA", (int(prefix_w + sub_w + suffix_w + 36 * render_scale), 70 * render_scale), (255, 255, 255, 0))
    ld = ImageDraw.Draw(label_img)
    ld.text((0, 0), prefix, fill=(25, 25, 25), font=y_axis_font)
    sub_x = prefix_w + 3
    ld.text((sub_x, 26 * render_scale), "2", fill=(25, 25, 25), font=sub_font)
    suffix_x = sub_x + sub_w + 8 * render_scale
    ld.text((suffix_x, 0), suffix, fill=(25, 25, 25), font=y_axis_font)
    rot = label_img.rotate(90, expand=True)
    img.paste(rot, (36 * render_scale, int(height / 2 - rot.height / 2)), rot)

    path.parent.mkdir(parents=True, exist_ok=True)
    if render_scale != 1:
        resampling = getattr(Image, "Resampling", Image).LANCZOS
        img = img.resize((final_width, final_height), resampling)
    img.save(path)


def alpha_grid(max_alpha: float = 256.0) -> np.ndarray:
    vals = np.unique(np.concatenate([
        np.linspace(1.05, 5.0, 12),
        np.linspace(5.5, 40.0, 28),
        np.linspace(45.0, max_alpha, 34),
    ]))
    return vals


def compute_full_identification(
    prior: AggregatePrior,
    sigmas: np.ndarray,
    ours_alphas: Tuple[float, ...],
    pac_grid: np.ndarray,
) -> Tuple[Dict[str, List[float]], List[dict]]:
    q_task = prior.prior_identification_success
    curves = {
        "Fano Mutual Information": [],
        "PAC-alpha(best)": [],
        **{f"PRW alpha={alpha:g}": [] for alpha in ours_alphas},
        "Ground truth": [],
    }
    rows = []
    for sigma in sigmas:
        mi = mutual_information(prior, float(sigma))
        pac_kl = math.log2(invert_kl_bound(q_task, mi))
        pac_alpha, best_alpha = pac_alpha_best_log2(prior, float(sigma), q_task, pac_grid)
        ours = []
        for alpha in ours_alphas:
            log_r = log_prw_alpha_cost(prior, float(sigma), alpha)
            ours.append(log_r / alpha / math.log(2.0))
        gt = log_ground_truth_full(prior, float(sigma)) / math.log(2.0)
        vals = [pac_kl, pac_alpha, *ours, gt]
        for key, val in zip(curves.keys(), vals):
            curves[key].append(val)
        row = {
            "figure": prior.name,
            "task": "full_identification",
            "sigma": sigma,
            "Fano_MI_log2_success": pac_kl,
            "PAC_alpha_best_log2_success": pac_alpha,
            "PAC_alpha_best_alpha": best_alpha,
            "ground_truth_log2_success": gt,
            "prior_success_log2": math.log2(q_task),
        }
        for alpha, value in zip(ours_alphas, ours):
            row[f"PRW_alpha_{alpha:g}_log2_success"] = value
        rows.append(row)
    return curves, rows


def compute_ball_task(
    prior: AggregatePrior,
    sigmas: np.ndarray,
    radius: int,
    ours_alphas: Tuple[float, ...],
    pac_grid: np.ndarray,
) -> Tuple[Dict[str, List[float]], List[dict]]:
    log_counts = ball_log_counts(prior.bits, radius)
    q_task = math.exp(log_prior_task_success(prior, radius))
    curves = {
        "Fano Mutual Information": [],
        "PAC-alpha(best)": [],
        **{f"PRW alpha={alpha:g}": [] for alpha in ours_alphas},
        "Ground truth": [],
    }
    rows = []
    for sigma in sigmas:
        mi = mutual_information(prior, float(sigma))
        pac_kl = math.log2(invert_kl_bound(q_task, mi))
        pac_alpha, best_alpha = pac_alpha_best_log2(prior, float(sigma), q_task, pac_grid)
        gt_log, prw_logs = ball_ground_truth_and_prw(prior, float(sigma), ours_alphas, log_counts)
        ours = [prw_logs[float(alpha)] / math.log(2.0) for alpha in ours_alphas]
        gt = gt_log / math.log(2.0)
        vals = [pac_kl, pac_alpha, *ours, gt]
        for key, val in zip(curves.keys(), vals):
            curves[key].append(val)
        row = {
            "figure": prior.name,
            "task": f"recover_at_least_{prior.bits - radius}_bits",
            "sigma": sigma,
            "Fano_MI_log2_success": pac_kl,
            "PAC_alpha_best_log2_success": pac_alpha,
            "PAC_alpha_best_alpha": best_alpha,
            "ground_truth_log2_success": gt,
            "prior_success_log2": math.log2(q_task),
        }
        for alpha, value in zip(ours_alphas, ours):
            row[f"PRW_alpha_{alpha:g}_log2_success"] = value
        rows.append(row)
    return curves, rows


def compute_zoom_full_identification(
    prior: AggregatePrior,
    sigmas: np.ndarray,
    ours_alphas: Tuple[float, ...],
    pac_grid: np.ndarray,
) -> Tuple[Dict[str, List[float]], List[dict]]:
    q_task = prior.prior_identification_success
    curves = {
        "PAC-alpha(best)": [],
        **{f"PRW alpha={alpha:g}": [] for alpha in ours_alphas},
    }
    rows = []
    for sigma in sigmas:
        pac_alpha, best_alpha = pac_alpha_best_log2(prior, float(sigma), q_task, pac_grid)
        ours = []
        for alpha in ours_alphas:
            log_r = log_prw_alpha_cost(prior, float(sigma), alpha)
            ours.append(log_r / alpha / math.log(2.0))
        vals = [pac_alpha, *ours]
        for key, val in zip(curves.keys(), vals):
            curves[key].append(val)
        row = {
            "figure": prior.name,
            "task": "full_identification_zoom",
            "sigma": sigma,
            "PAC_alpha_best_log2_success": pac_alpha,
            "PAC_alpha_best_alpha": best_alpha,
            "prior_success_log2": math.log2(q_task),
        }
        for alpha, value in zip(ours_alphas, ours):
            row[f"PRW_alpha_{alpha:g}_log2_success"] = value
        rows.append(row)
    return curves, rows


def compute_zoom_ball_task(
    prior: AggregatePrior,
    sigmas: np.ndarray,
    radius: int,
    ours_alphas: Tuple[float, ...],
    pac_grid: np.ndarray,
) -> Tuple[Dict[str, List[float]], List[dict]]:
    log_counts = ball_log_counts(prior.bits, radius)
    q_task = math.exp(log_prior_task_success(prior, radius))
    curves = {
        "PAC-alpha(best)": [],
        **{f"PRW alpha={alpha:g}": [] for alpha in ours_alphas},
    }
    rows = []
    for sigma in sigmas:
        pac_alpha, best_alpha = pac_alpha_best_log2(prior, float(sigma), q_task, pac_grid)
        prw_logs = ball_prw_only(prior, float(sigma), ours_alphas, log_counts)
        ours = [prw_logs[float(alpha)] / math.log(2.0) for alpha in ours_alphas]
        vals = [pac_alpha, *ours]
        for key, val in zip(curves.keys(), vals):
            curves[key].append(val)
        row = {
            "figure": prior.name,
            "task": f"recover_at_least_{prior.bits - radius}_bits_zoom",
            "sigma": sigma,
            "PAC_alpha_best_log2_success": pac_alpha,
            "PAC_alpha_best_alpha": best_alpha,
            "prior_success_log2": math.log2(q_task),
        }
        for alpha, value in zip(ours_alphas, ours):
            row[f"PRW_alpha_{alpha:g}_log2_success"] = value
        rows.append(row)
    return curves, rows


def csv_task_name(spec: dict) -> str:
    prior = spec["prior"]
    if spec["task"] == "full":
        return "full_identification"
    return f"recover_at_least_{prior.bits - spec['radius']}_bits"


def load_cached_curves(
    cache_path: Path,
    spec: dict,
    sigmas: np.ndarray,
    ours_alphas: Tuple[float, ...],
) -> Optional[Tuple[Dict[str, List[float]], List[dict]]]:
    if not cache_path.exists():
        return None
    figure_name = spec["prior"].name
    task_name = csv_task_name(spec)
    with cache_path.open() as f:
        rows = [
            row for row in csv.DictReader(f)
            if row.get("figure") == figure_name and row.get("task") == task_name
        ]
    if len(rows) != len(sigmas):
        return None
    rows.sort(key=lambda row: float(row["sigma"]))
    cached_sigmas = np.asarray([float(row["sigma"]) for row in rows])
    if not np.allclose(cached_sigmas, sigmas):
        return None

    curves = {
        "Fano Mutual Information": [float(row["Fano_MI_log2_success"]) for row in rows],
        "PAC-alpha(best)": [float(row["PAC_alpha_best_log2_success"]) for row in rows],
    }
    for alpha in ours_alphas:
        curves[f"PRW alpha={alpha:g}"] = [float(row[f"PRW_alpha_{alpha:g}_log2_success"]) for row in rows]
    curves["Ground truth"] = [float(row["ground_truth_log2_success"]) for row in rows]
    return curves, rows


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows: List[dict] = []
    zoom_rows: List[dict] = []
    sigmas_025_to_30 = np.arange(0.25, 30.0 + 1e-9, 0.25)
    sigmas_40_to_80 = np.arange(40.0, 80.0 + 1e-9, 1.0)
    prw_alphas = (16.0, 64.0, 256.0)
    pac_grid = np.unique(np.concatenate([
        np.linspace(1.001, 10.0, 181),
        np.linspace(10.25, 64.0, 108),
        np.linspace(68.0, 256.0, 48),
    ]))

    specs = [
        {
            "filename": "fig1_uniform128_full_sigma025_30.png",
            "title": "Uniform 128-Bit Secret: Full Identification",
            "panel_title": "(a) Uniform 128-bit, full identification",
            "subtitle": "",
            "prior": uniform_prior(128),
            "sigmas": sigmas_025_to_30,
            "ours_alphas": prw_alphas,
            "task": "full",
        },
        {
            "filename": "fig2_uniform256_full_sigma025_30.png",
            "title": "Uniform 256-Bit Secret: Full Identification",
            "panel_title": "(b) Uniform 256-bit, full identification",
            "subtitle": "",
            "prior": uniform_prior(256),
            "sigmas": sigmas_025_to_30,
            "ours_alphas": prw_alphas,
            "task": "full",
        },
        {
            "filename": "fig3_poisson_hamming_weight96_256_full_sigma025_30.png",
            "title": "Poisson-HW(λ=96) 256-Bit Secret: Full Identification",
            "panel_title": "(c) Poisson-HW 256-bit full identification",
            "subtitle": "",
            "prior": hamming_poisson_prior(256, lam=96.0),
            "sigmas": sigmas_025_to_30,
            "ours_alphas": prw_alphas,
            "task": "full",
        },
        {
            "filename": "fig4_uniform256_recover240_sigma025_30.png",
            "title": "Uniform 256-Bit Secret: Recover at Least 240 Bits",
            "panel_title": "(d) Uniform 256-bit, recover at least 240 bit",
            "subtitle": "",
            "prior": uniform_prior(256),
            "sigmas": sigmas_025_to_30,
            "ours_alphas": prw_alphas,
            "task": "ball",
            "radius": 16,
        },
    ]

    combined_panels = []
    cache_path = OUT / "gaussian_comparison_results.csv"
    for spec in specs:
        prior = spec["prior"]
        cached = load_cached_curves(cache_path, spec, spec["sigmas"], spec["ours_alphas"])
        if cached is not None:
            curves, rows = cached
            print(f"loaded cached main curves for {spec['filename']}", flush=True)
        elif spec["task"] == "full":
            curves, rows = compute_full_identification(prior, spec["sigmas"], spec["ours_alphas"], pac_grid)
        else:
            curves, rows = compute_ball_task(prior, spec["sigmas"], spec["radius"], spec["ours_alphas"], pac_grid)

        if spec["task"] == "full":
            zoom_curves, zoom_rows_for_spec = compute_zoom_full_identification(prior, sigmas_40_to_80, spec["ours_alphas"], pac_grid)
        else:
            zoom_curves, zoom_rows_for_spec = compute_zoom_ball_task(prior, sigmas_40_to_80, spec["radius"], spec["ours_alphas"], pac_grid)
        all_rows.extend(rows)
        zoom_rows.extend(zoom_rows_for_spec)
        combined_panels.append((spec["panel_title"], spec["sigmas"], curves, sigmas_40_to_80, zoom_curves))
        plot_curves(OUT / spec["filename"], spec["title"], spec["subtitle"], spec["sigmas"], curves)
        print(f"wrote {spec['filename']}", flush=True)
    plot_combined_curves(OUT / "fig1_4_combined_1x4_sigma025_30_zoom4080_v7.png", combined_panels)
    print("wrote fig1_4_combined_1x4_sigma025_30_zoom4080_v7.png", flush=True)

    fieldnames: List[str] = []
    for row in all_rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with (OUT / "gaussian_comparison_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    zoom_fieldnames: List[str] = []
    for row in zoom_rows:
        for key in row:
            if key not in zoom_fieldnames:
                zoom_fieldnames.append(key)
    with (OUT / "gaussian_comparison_zoom_40_80_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=zoom_fieldnames)
        writer.writeheader()
        writer.writerows(zoom_rows)
    print(f"wrote results to {OUT}")


if __name__ == "__main__":
    run()
