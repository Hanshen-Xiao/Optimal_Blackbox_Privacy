"""Discrete-noise experiments for black-box privatization of entropic secrets.

The experiments use finite output domains and discrete additive noise.  The
large secret spaces are handled by product-channel compression:

* AES power: one coordinate is one secret byte.  The leakage is the Hamming
  weight of an AES S-box output, so the per-coordinate leakage alphabet is
  {0,...,8}.  A secret with l bits has l/8 independent coordinates.

* RSA timing: one coordinate is one exponent bit observed through a discretized
  square-and-multiply timing window.  The leakage alphabet is {0,1}; a noise bin
  is interpreted as 1 ns by default.

For each coordinate channel we calibrate a discrete additive-noise family by
binary search over a scale parameter, minimizing the second moment within that
family.  The finite support is explicit and all accounting quantities are exact
on the compressed channel.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont


Array = np.ndarray
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "discrete"


@dataclass
class Channel:
    name: str
    leakage_values: Array
    leakage_probs: Array
    sum_prior_alpha: Callable[[float], Array]
    max_prior_by_leakage: Array
    prior_success_coord: float
    coord_bits: int
    unit_label: str
    unit_scale: float = 1.0


@dataclass
class Calibration:
    status: str
    scale: float
    second_moment: float
    mean: float
    log_bound_total: float
    support_size: int


def n_choose_k_counts(bits: int) -> Array:
    return np.array([math.comb(bits, k) for k in range(bits + 1)], dtype=float)


def make_aes_power_channel() -> Channel:
    counts = n_choose_k_counts(8)
    probs = counts / 256.0
    values = np.arange(9, dtype=int)

    def sum_prior_alpha(alpha: float) -> Array:
        return counts * (1.0 / 256.0) ** alpha

    return Channel(
        name="AES power",
        leakage_values=values,
        leakage_probs=probs,
        sum_prior_alpha=sum_prior_alpha,
        max_prior_by_leakage=np.where(counts > 0, 1.0 / 256.0, 0.0),
        prior_success_coord=1.0 / 256.0,
        coord_bits=8,
        unit_label="power bin",
    )


def make_rsa_timing_bit_channel(bit_prob_one: float = 0.5, unit_ns: float = 1.0) -> Channel:
    p1 = float(bit_prob_one)
    probs_by_secret = np.array([1.0 - p1, p1], dtype=float)
    values = np.array([0, 1], dtype=int)
    leakage_probs = probs_by_secret.copy()

    def sum_prior_alpha(alpha: float) -> Array:
        return probs_by_secret**alpha

    return Channel(
        name="RSA timing",
        leakage_values=values,
        leakage_probs=leakage_probs,
        sum_prior_alpha=sum_prior_alpha,
        max_prior_by_leakage=probs_by_secret.copy(),
        prior_success_coord=float(np.max(probs_by_secret)),
        coord_bits=1,
        unit_label="ns",
        unit_scale=unit_ns,
    )


def support(kind: str, k: int) -> Array:
    if kind == "symmetric":
        return np.arange(-k, k + 1, dtype=int)
    if kind == "onesided":
        return np.arange(0, k + 1, dtype=int)
    raise ValueError(f"unknown noise kind: {kind}")


def discrete_gaussian_noise(kind: str, k: int, scale: float) -> Tuple[Array, Array]:
    e = support(kind, k)
    if scale <= 0:
        q = np.zeros_like(e, dtype=float)
        q[np.where(e == 0)[0][0]] = 1.0
        return e, q
    logq = -(e.astype(float) ** 2) / (2.0 * scale * scale)
    logq -= float(np.max(logq))
    q = np.exp(logq)
    q /= float(np.sum(q))
    return e, q


def second_moment(e: Array, q: Array) -> float:
    return float(np.sum(q * (e.astype(float) ** 2)))


def mean_noise(e: Array, q: Array) -> float:
    return float(np.sum(q * e.astype(float)))


def shifted_noise_prob(y: int, h: int, e_to_q: Dict[int, float]) -> float:
    return e_to_q.get(int(y - h), 0.0)


def y_domain(channel: Channel, e: Array) -> range:
    return range(int(channel.leakage_values[0] + e[0]), int(channel.leakage_values[-1] + e[-1]) + 1)


def log_prw_coord(channel: Channel, e: Array, q: Array, alpha: float) -> float:
    e_to_logq = {int(v): math.log(float(p)) for v, p in zip(e, q) if p > 0.0}
    s_alpha = channel.sum_prior_alpha(alpha)
    log_s_alpha = np.full_like(s_alpha, float("-inf"), dtype=float)
    positive = s_alpha > 0
    log_s_alpha[positive] = np.log(s_alpha[positive])
    log_terms = []
    for y in y_domain(channel, e):
        inner_logs = []
        for idx, h in enumerate(channel.leakage_values):
            logq = e_to_logq.get(int(y - h))
            if logq is not None and math.isfinite(float(log_s_alpha[idx])):
                inner_logs.append(float(log_s_alpha[idx]) + alpha * logq)
        if inner_logs:
            m = max(inner_logs)
            log_inner = m + math.log(sum(math.exp(v - m) for v in inner_logs))
            log_terms.append(log_inner / alpha)
    if not log_terms:
        return float("-inf")
    m = max(log_terms)
    return m + math.log(sum(math.exp(v - m) for v in log_terms))


def log_exact_coord(channel: Channel, e: Array, q: Array) -> float:
    e_to_q = {int(v): float(p) for v, p in zip(e, q) if p > 0.0}
    total = 0.0
    for y in y_domain(channel, e):
        best = 0.0
        for idx, h in enumerate(channel.leakage_values):
            p = shifted_noise_prob(y, int(h), e_to_q)
            best = max(best, float(channel.max_prior_by_leakage[idx]) * p)
        total += best
    return math.log(total) if total > 0 else float("-inf")


def entropy(probs: Array) -> float:
    p = probs[probs > 0]
    return float(-np.sum(p * np.log(p)))


def mutual_information_coord(channel: Channel, e: Array, q: Array) -> float:
    e_to_q = {int(v): float(p) for v, p in zip(e, q) if p > 0.0}
    py = []
    for y in y_domain(channel, e):
        p = 0.0
        for h, ph in zip(channel.leakage_values, channel.leakage_probs):
            p += float(ph) * shifted_noise_prob(y, int(h), e_to_q)
        py.append(p)
    return entropy(np.array(py, dtype=float)) - entropy(q)


def binary_kl(p: float, q: float) -> float:
    eps = 1e-300
    p = min(max(p, eps), 1.0 - eps)
    q = min(max(q, eps), 1.0 - eps)
    return p * math.log(p / q) + (1.0 - p) * math.log((1.0 - p) / (1.0 - q))


def full_log_bound(
    channel: Channel,
    coord_count: int,
    e: Array,
    q: Array,
    method: str,
    alpha: Optional[float] = None,
) -> float:
    if method == "exact":
        return coord_count * log_exact_coord(channel, e, q)
    if method == "prw":
        if alpha is None:
            raise ValueError("alpha is required for PRW")
        return coord_count * log_prw_coord(channel, e, q, alpha)
    raise ValueError(method)


def calibrate_noise(
    channel: Channel,
    coord_count: int,
    target_log2: float,
    kind: str,
    method: str,
    alpha: Optional[float],
    k: int,
    scale_hi: float,
) -> Calibration:
    target_log = target_log2 * math.log(2.0)
    e0, q0 = discrete_gaussian_noise(kind, k, 0.0)

    if method in {"exact", "prw"}:
        raw = full_log_bound(channel, coord_count, e0, q0, method, alpha)
        if raw <= target_log:
            return Calibration("ok", 0.0, 0.0, 0.0, raw, len(e0))

        def predicate(scale: float) -> Tuple[bool, float, Array, Array]:
            e, q = discrete_gaussian_noise(kind, k, scale)
            val = full_log_bound(channel, coord_count, e, q, method, alpha)
            return val <= target_log, val, e, q

    elif method == "kl":
        prior_success = channel.prior_success_coord**coord_count
        target = 2.0**target_log2
        if target <= prior_success:
            return Calibration("target_below_prior", 0.0, 0.0, 0.0, math.log(prior_success), len(e0))
        kl_threshold = binary_kl(target, prior_success)

        raw_mi = coord_count * mutual_information_coord(channel, e0, q0)
        if raw_mi <= kl_threshold:
            return Calibration("ok", 0.0, 0.0, 0.0, raw_mi, len(e0))

        def predicate(scale: float) -> Tuple[bool, float, Array, Array]:
            e, q = discrete_gaussian_noise(kind, k, scale)
            val = coord_count * mutual_information_coord(channel, e, q)
            return val <= kl_threshold, val, e, q

    else:
        raise ValueError(method)

    ok, val, e, q = predicate(scale_hi)
    if not ok:
        return Calibration("not_reached", scale_hi, second_moment(e, q), mean_noise(e, q), val, len(e))

    lo, hi = 0.0, scale_hi
    best_val, best_e, best_q = val, e, q
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        ok, val, e, q = predicate(mid)
        if ok:
            hi = mid
            best_val, best_e, best_q = val, e, q
        else:
            lo = mid
    return Calibration(
        "ok",
        hi,
        second_moment(best_e, best_q),
        mean_noise(best_e, best_q),
        best_val,
        len(best_e),
    )


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
    log_x: bool = False,
    log_y: bool = False,
) -> None:
    w, h = 1450, 900
    m = {"left": 145, "right": 260, "top": 115, "bottom": 120}
    pw, ph = w - m["left"] - m["right"], h - m["top"] - m["bottom"]
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    ft, fs, fa, fk = font(32, True), font(18), font(23), font(18)
    d.text((m["left"], 32), title, fill=(25, 25, 25), font=ft)
    d.text((m["left"], 75), subtitle, fill=(80, 80, 80), font=fs)
    all_x = [x for pts in series.values() for x, _ in pts]
    all_y = [y for pts in series.values() for _, y in pts if y > 0]
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if log_y:
        ymin = 10 ** math.floor(math.log10(max(ymin, 1e-12)))
        ymax = 10 ** math.ceil(math.log10(ymax))
    else:
        ymin = 0.0
        ymax *= 1.12

    def xm(x: float) -> float:
        if log_x:
            return m["left"] + (math.log10(max(x, xmin)) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin)) * pw
        return m["left"] + (x - xmin) / max(1e-12, xmax - xmin) * pw

    def ym(y: float) -> float:
        if log_y:
            return m["top"] + (math.log10(ymax) - math.log10(max(y, ymin))) / (math.log10(ymax) - math.log10(ymin)) * ph
        return m["top"] + (ymax - y) / (ymax - ymin) * ph

    d.rectangle([m["left"], m["top"], m["left"] + pw, m["top"] + ph], fill=(250, 250, 250), outline=(180, 180, 180))
    if log_x:
        x_ticks = []
        for exp in range(int(math.log10(xmin)), int(math.log10(xmax)) + 1):
            for mult in [1, 2, 5]:
                val = mult * 10**exp
                if xmin <= val <= xmax:
                    x_ticks.append(val)
    else:
        x_ticks = [xmin + (xmax - xmin) * i / 5 for i in range(6)]
    for x in x_ticks:
        px = xm(x)
        d.line([(px, m["top"]), (px, m["top"] + ph)], fill=(225, 225, 225))
        label = f"{x:g}"
        tw = d.textlength(label, font=fk)
        d.text((px - tw / 2, m["top"] + ph + 14), label, fill=(45, 45, 45), font=fk)
    if log_y:
        y_ticks = []
        e0, e1 = int(math.log10(ymin)), int(math.log10(ymax))
        for exp in range(e0, e1 + 1):
            for mult in [1, 2, 5]:
                val = mult * 10**exp
                if ymin <= val <= ymax:
                    y_ticks.append(val)
    else:
        y_ticks = [ymin + (ymax - ymin) * i / 5 for i in range(6)]
    for y in y_ticks:
        py = ym(y)
        d.line([(m["left"], py), (m["left"] + pw, py)], fill=(225, 225, 225))
        label = f"{y:.2g}" if y < 1000 else f"{y/1000:.1f}k"
        tw = d.textlength(label, font=fk)
        d.text((m["left"] - tw - 12, py - 9), label, fill=(45, 45, 45), font=fk)
    d.line([(m["left"], m["top"] + ph), (m["left"] + pw, m["top"] + ph)], fill=(30, 30, 30), width=2)
    d.line([(m["left"], m["top"]), (m["left"], m["top"] + ph)], fill=(30, 30, 30), width=2)

    colors = [(0, 114, 178), (213, 94, 0), (0, 158, 115), (204, 121, 167), (86, 180, 233), (230, 159, 0)]
    for idx, (name, pts) in enumerate(series.items()):
        color = colors[idx % len(colors)]
        xy = [(xm(x), ym(max(y, ymin))) for x, y in pts if y > 0]
        for a, b in zip(xy, xy[1:]):
            d.line([a, b], fill=color, width=4)
        for px, py in xy:
            d.ellipse([px - 6, py - 6, px + 6, py + 6], fill=color, outline="white", width=2)
        ly = m["top"] + 20 + idx * 32
        lx = m["left"] + pw + 35
        d.line([(lx, ly + 10), (lx + 45, ly + 10)], fill=color, width=4)
        d.text((lx + 58, ly), name, fill=(35, 35, 35), font=fk)

    tw = d.textlength(xlabel, font=fa)
    d.text((m["left"] + pw / 2 - tw / 2, h - 70), xlabel, fill=(30, 30, 30), font=fa)
    label_img = Image.new("RGBA", (700, 38), (255, 255, 255, 0))
    ld = ImageDraw.Draw(label_img)
    ld.text((0, 0), ylabel, fill=(30, 30, 30), font=fa)
    rot = label_img.rotate(90, expand=True)
    img.paste(rot, (30, m["top"] + ph // 2 - rot.height // 2), rot)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    aes = make_aes_power_channel()
    rsa = make_rsa_timing_bit_channel(unit_ns=1.0)
    rows: List[dict] = []
    k = 512
    scale_hi = 600.0

    # 1. AES entropy and target sweep.
    entropies = [208, 224, 240, 256]
    targets = [-200, -150, -100]
    for l in entropies:
        n = l // aes.coord_bits
        for target in targets:
            cal = calibrate_noise(aes, n, target, "symmetric", "prw", 128.0, k, scale_hi)
            rows.append({
                "experiment": "aes_entropy",
                "channel": aes.name,
                "entropy_bits": l,
                "target_log2": target,
                "constraint": "zero_mean",
                "method": "PRW alpha=128",
                "status": cal.status,
                "second_moment": cal.second_moment,
                "mean_noise": cal.mean,
                "scale": cal.scale,
            })

    # 2. RSA target/constraint sweep.
    for target in targets:
        for kind, label in [("symmetric", "zero_mean"), ("onesided", "delay_only")]:
            cal = calibrate_noise(rsa, 256, target, kind, "prw", 512.0, k, scale_hi)
            rows.append({
                "experiment": "rsa_constraint",
                "channel": rsa.name,
                "entropy_bits": 256,
                "target_log2": target,
                "constraint": label,
                "method": "PRW alpha=512",
                "status": cal.status,
                "second_moment": cal.second_moment,
                "rms_ns_per_bit": math.sqrt(cal.second_moment),
                "mean_delay_ns_per_bit": cal.mean,
                "mean_delay_ns_total": cal.mean * 256,
                "scale": cal.scale,
            })

    # 3. Alpha convergence against exact ground truth for RSA.
    for alpha in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        cal = calibrate_noise(rsa, 256, -150, "symmetric", "prw", float(alpha), k, scale_hi)
        rows.append({
            "experiment": "alpha_convergence",
            "channel": rsa.name,
            "entropy_bits": 256,
            "target_log2": -150,
            "alpha": alpha,
            "method": f"PRW alpha={alpha}",
            "constraint": "zero_mean",
            "status": cal.status,
            "second_moment": cal.second_moment,
            "scale": cal.scale,
        })
    exact_cal = calibrate_noise(rsa, 256, -150, "symmetric", "exact", None, k, scale_hi)
    rows.append({
        "experiment": "alpha_convergence",
        "channel": rsa.name,
        "entropy_bits": 256,
        "target_log2": -150,
        "alpha": "exact",
        "method": "Exact MAP",
        "constraint": "zero_mean",
        "status": exact_cal.status,
        "second_moment": exact_cal.second_moment,
        "scale": exact_cal.scale,
    })

    # 4. KL/PAC comparison, uniform and non-uniform RSA bit secret.
    rsa_biased = make_rsa_timing_bit_channel(bit_prob_one=0.6, unit_ns=1.0)
    for channel, label, target in [(rsa, "uniform", -150), (rsa_biased, "nonuniform_p=0.6", -120)]:
        for method, alpha, method_label in [("exact", None, "Exact MAP"), ("prw", 512.0, "PRW alpha=512"), ("kl", None, "PAC KL/MI")]:
            cal = calibrate_noise(channel, 256, target, "symmetric", method, alpha, k, scale_hi)
            rows.append({
                "experiment": "kl_comparison",
                "prior": label,
                "channel": channel.name,
                "entropy_bits": 256 if label == "uniform" else "biased product",
                "target_log2": target,
                "constraint": "zero_mean",
                "method": method_label,
                "status": cal.status,
                "second_moment": cal.second_moment,
                "scale": cal.scale,
            })

    # 5. Composition: fixed total posterior success, equal budget split.
    l, alpha, target = 256, 512.0, -150
    coord_count = l
    base_coord = rsa.sum_prior_alpha(alpha).sum()
    log_base_full = coord_count * math.log(base_coord)
    log_target_alpha = alpha * target * math.log(2.0)
    total_log_r = log_target_alpha - log_base_full
    for T in range(1, 9):
        log_r_round_full = total_log_r / T
        log_r_byte = log_r_round_full / coord_count
        # Need PRW coord delta satisfying delta_coord^alpha/base_coord <= exp(log_r_byte).
        coord_target_log = (math.log(base_coord) + log_r_byte) / alpha
        pseudo_target_log2_full = coord_target_log / math.log(2.0) * coord_count
        cal = calibrate_noise(rsa, coord_count, pseudo_target_log2_full, "onesided", "prw", alpha, k, scale_hi)
        rows.append({
            "experiment": "composition",
            "channel": rsa.name,
            "entropy_bits": l,
            "target_log2": target,
            "constraint": "delay_only",
            "method": "PRW alpha=512",
            "T": T,
            "status": cal.status,
            "avg_second_moment_per_iter": cal.second_moment,
            "mean_delay_ns_per_iter": cal.mean * l,
            "scale": cal.scale,
        })

    write_csv(OUT / "discrete_experiment_results.csv", rows)

    # Plots.
    aes_series: Dict[str, List[Tuple[float, float]]] = {}
    for target in targets:
        pts = []
        for r in rows:
            if r.get("experiment") == "aes_entropy" and r.get("target_log2") == target:
                pts.append((float(r["entropy_bits"]), max(float(r["second_moment"]), 1e-9)))
        aes_series[f"target 2^{target}"] = sorted(pts)
    line_plot(
        OUT / "aes_entropy_targets.png",
        "AES Power: Noise vs. Secret Entropy",
        "Discrete zero-mean noise, PRW alpha=128, objective E[N^2]",
        "Secret entropy l (bits)",
        "Noise second moment",
        aes_series,
        log_x=False,
        log_y=True,
    )

    rsa_series: Dict[str, List[Tuple[float, float]]] = {"zero-mean": [], "delay-only": []}
    for r in rows:
        if r.get("experiment") == "rsa_constraint":
            key = "zero-mean" if r["constraint"] == "zero_mean" else "delay-only"
            rsa_series[key].append((-float(r["target_log2"]), max(float(r["second_moment"]), 1e-9)))
    line_plot(
        OUT / "rsa_noise_constraints.png",
        "RSA Timing: Constraint Cost",
        "Discrete timing noise, PRW alpha=512; delay-only overhead is nonnegative ns bins",
        "Required posterior exponent s for target 2^-s",
        "Per-bit noise second moment (ns^2)",
        {k: sorted(v) for k, v in rsa_series.items()},
        log_x=False,
        log_y=True,
    )

    alpha_pts = []
    exact_y = max(exact_cal.second_moment, 1e-9)
    exact_series = []
    for r in rows:
        if (
            r.get("experiment") == "alpha_convergence"
            and isinstance(r.get("alpha"), int)
            and r.get("status") == "ok"
        ):
            alpha_pts.append((float(r["alpha"]), max(float(r["second_moment"]), 1e-9)))
            exact_series.append((float(r["alpha"]), exact_y))
    line_plot(
        OUT / "rsa_alpha_convergence.png",
        "RSA Timing: Alpha Convergence",
        "Target 2^-150, zero-mean discrete noise; exact MAP shown as reference",
        "alpha",
        "Per-bit noise second moment (ns^2)",
        {"PRW": sorted(alpha_pts), "Exact MAP": sorted(exact_series)},
        log_x=True,
        log_y=True,
    )

    comp_pts = []
    for r in rows:
        if r.get("experiment") == "composition":
            comp_pts.append((float(r["T"]), max(float(r["avg_second_moment_per_iter"]), 1e-9)))
    line_plot(
        OUT / "rsa_composition_delay.png",
        "RSA Timing Composition",
        "Fixed total target 2^-150, one-sided delay noise, PRW alpha=512",
        "Composition rounds T",
        "Average per-round second moment (ns^2)",
        {"delay-only": sorted(comp_pts)},
        log_x=False,
        log_y=True,
    )

    # KL comparison as a compact bar-like line over methods.
    kl_series: Dict[str, List[Tuple[float, float]]] = {}
    method_order = {"Exact MAP": 1, "PRW alpha=512": 2, "PAC KL/MI": 3}
    for prior in ["uniform", "nonuniform_p=0.6"]:
        pts = []
        for r in rows:
            if r.get("experiment") == "kl_comparison" and r.get("prior") == prior:
                y = float(r["second_moment"])
                if r["status"] == "not_reached":
                    y = max(y, 1e5)
                pts.append((method_order[str(r["method"])], max(y, 1e-9)))
        kl_series[prior] = sorted(pts)
    line_plot(
        OUT / "rsa_kl_comparison.png",
        "PAC KL/MI Comparison",
        "PAC KL bars at cap indicate target not reached within support [-512,512]",
        "Method index: 1 Exact, 2 PRW, 3 PAC KL/MI",
        "Per-bit noise second moment (ns^2)",
        kl_series,
        log_x=False,
        log_y=True,
    )

    print(f"wrote results and plots to {OUT}")


if __name__ == "__main__":
    run()
