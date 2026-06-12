"""Render the eight experiment panels as one 2x4 figure.

This script reads the already-computed CSV outputs from the two experiment
generators and draws a single publication-style panel figure.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from PIL import Image, ImageDraw

from gaussian_comparison_figures import display_curve_label, font


ROOT = Path(__file__).resolve().parent
GAUSSIAN_CSV = ROOT / "outputs" / "gaussian_comparison" / "gaussian_comparison_results.csv"
GAUSSIAN_ZOOM_CSV = ROOT / "outputs" / "gaussian_comparison" / "gaussian_comparison_zoom_40_80_results.csv"
MECHANISM_CSV = ROOT / "outputs" / "figures_5_to_8" / "figures_5_to_8_results.csv"
OUT = ROOT / "outputs" / "combined"


PointSeries = Dict[str, List[Tuple[float, float]]]


@dataclass
class Panel:
    title: str
    series: PointSeries
    x_ticks: Sequence[float]
    y_label: str = ""
    y_zero: bool = False
    draw_markers: bool = False
    inline_legend: bool = False
    color_group: str = "gaussian"
    zoom_series: Optional[PointSeries] = None


GAUSSIAN_ORDER = [
    "Fano Mutual Information",
    "PAC-alpha(best)",
    "PRW alpha=16",
    "PRW alpha=64",
    "PRW alpha=256",
    "Ground truth",
]

MECHANISM_ORDER = [
    "Gaussian",
    "Optimal zero-mean",
    "Optimal non-negative noise",
]

GAUSSIAN_COLORS = {
    "Fano Mutual Information": (0, 114, 178),
    "PAC-alpha(best)": (213, 94, 0),
    "PRW alpha=16": (0, 158, 115),
    "PRW alpha=64": (230, 159, 0),
    "PRW alpha=256": (204, 121, 167),
    "Ground truth": (150, 150, 150),
}

ENTROPY_COLORS = {
    "l=253": (0, 114, 178),
    "l=254": (213, 94, 0),
    "l=255": (0, 158, 115),
    "l=256": (204, 121, 167),
}

MECHANISM_COLORS = {
    "Gaussian": (86, 180, 233),
    "Optimal zero-mean": (117, 112, 179),
    "Optimal non-negative noise": (166, 86, 40),
}


def load_rows(path: Path) -> List[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def as_float(value: str) -> float | None:
    if value == "":
        return None
    return float(value)


def gaussian_panels(rows: Sequence[dict], zoom_rows: Sequence[dict]) -> List[Panel]:
    specs = [
        ("(a) Uniform 128-bit, full identification", "uniform 128-bit", "full_identification"),
        ("(b) Uniform 256-bit, full identification", "uniform 256-bit", "full_identification"),
        ("(c) Poisson-HW 256-bit full identification", "Poisson-HW lambda=96", "full_identification"),
        ("(d) Uniform 256-bit, recover at least 240 bit", "uniform 256-bit", "recover_at_least_240_bits"),
    ]
    field_for_name = {
        "Fano Mutual Information": "Fano_MI_log2_success",
        "PAC-alpha(best)": "PAC_alpha_best_log2_success",
        "PRW alpha=16": "PRW_alpha_16_log2_success",
        "PRW alpha=64": "PRW_alpha_64_log2_success",
        "PRW alpha=256": "PRW_alpha_256_log2_success",
        "Ground truth": "ground_truth_log2_success",
    }
    zoom_field_for_name = {
        "PAC-alpha(best)": "PAC_alpha_best_log2_success",
        "PRW alpha=16": "PRW_alpha_16_log2_success",
        "PRW alpha=64": "PRW_alpha_64_log2_success",
        "PRW alpha=256": "PRW_alpha_256_log2_success",
    }
    panels: List[Panel] = []
    for title, figure_name, task in specs:
        selected = [row for row in rows if row["figure"] == figure_name and row["task"] == task]
        selected.sort(key=lambda row: float(row["sigma"]))
        series: PointSeries = {}
        for name in GAUSSIAN_ORDER:
            field = field_for_name[name]
            series[name] = [(float(row["sigma"]), float(row[field])) for row in selected]
        zoom_task = f"{task}_zoom"
        selected_zoom = [
            row for row in zoom_rows
            if row["figure"] == figure_name and row["task"] == zoom_task
        ]
        selected_zoom.sort(key=lambda row: float(row["sigma"]))
        zoom_series: PointSeries = {}
        for name, field in zoom_field_for_name.items():
            zoom_series[name] = [(float(row["sigma"]), float(row[field])) for row in selected_zoom]
        panels.append(
            Panel(
                title=title,
                series=series,
                x_ticks=[0.25, 10, 20, 30],
                color_group="gaussian",
                zoom_series=zoom_series,
            )
        )
    return panels


def mechanism_panels(rows: Sequence[dict]) -> List[Panel]:
    series5: PointSeries = {}
    for bits in [253, 254, 255, 256]:
        selected = [
            row for row in rows
            if row["figure"] == "5" and row["mechanism"] == "Gaussian" and row["secret_bits"] == str(bits)
        ]
        selected.sort(key=lambda row: float(row["target_exponent"]))
        series5[f"l={bits}"] = [
            (float(row["target_exponent"]), float(row["y_value"])) for row in selected
        ]

    panels = [
        Panel(
            title="(e) Effect from various $l$-bit secret entropy",
            series=series5,
            x_ticks=[248, 249, 250, 251, 252],
            y_label="Gaussian noise std σ",
            y_zero=True,
            draw_markers=False,
            inline_legend=True,
            color_group="entropy",
        )
    ]

    series6: PointSeries = {name: [] for name in MECHANISM_ORDER}
    mechanism_map = {
        "Gaussian": "Gaussian",
        "optimal_discrete_symmetric": "Optimal zero-mean",
        "optimal_discrete_nonnegative": "Optimal non-negative noise",
    }
    selected6 = [row for row in rows if row["figure"] == "6"]
    for row in selected6:
        label = mechanism_map[row["mechanism"]]
        y_value = as_float(row["expected_abs_noise"])
        if y_value is not None:
            series6[label].append((float(row["target_exponent"]), y_value))
    for values in series6.values():
        values.sort()
    panels.append(
        Panel(
            title="(f) Randomization mechanisms",
            series=series6,
            x_ticks=[248, 249, 250, 251, 252, 253, 254],
            y_label="Expected absolute perturbation",
            y_zero=True,
            draw_markers=False,
            color_group="mechanism",
        )
    )

    for fig_no, exponent in [(7, 240), (8, 220)]:
        series: PointSeries = {name: [] for name in MECHANISM_ORDER}
        selected = [row for row in rows if row["figure"] == str(fig_no)]
        for row in selected:
            label = mechanism_map[row["mechanism"]]
            y_value = as_float(row["expected_abs_noise_per_round"])
            if y_value is not None:
                series[label].append((float(row["rounds"]), y_value))
        for values in series.values():
            values.sort()
        panels.append(
            Panel(
                title=f"({chr(ord('a') + fig_no - 1)}) Composition, full identification $2^{{-{exponent}}}$",
                series=series,
                x_ticks=[1, 10, 20, 30, 40, 50],
                y_label="Average perturbation up to iteration t" if fig_no == 7 else "",
                y_zero=True,
                draw_markers=True,
                color_group="mechanism",
            )
        )
    return panels


def color_for(panel: Panel, name: str) -> Tuple[int, int, int]:
    if panel.color_group == "entropy":
        return ENTROPY_COLORS[name]
    if panel.color_group == "mechanism":
        return MECHANISM_COLORS[name]
    return GAUSSIAN_COLORS[name]


def draw_rich_text(
    d: ImageDraw.ImageDraw,
    xy: Tuple[float, float],
    text: str,
    base_font,
    *,
    scale: int,
    fill: Tuple[int, int, int] = (25, 25, 25),
) -> float:
    """Draw a small subset of mathtext used by the figure titles."""
    x, y = xy
    small_font = font(max(1, int(base_font.size * 0.78)), True)
    i = 0
    while i < len(text):
        if text.startswith("$l$", i):
            d.text((x, y), "l", fill=fill, font=base_font)
            x += d.textlength("l", font=base_font)
            i += 3
            continue
        if text.startswith("$2^{-", i):
            end = text.find("}$", i)
            if end != -1:
                exponent = "-" + text[i + 5:end]
                d.text((x, y), "2", fill=fill, font=base_font)
                x += d.textlength("2", font=base_font) + 1 * scale
                d.text((x, y - 12 * scale), exponent, fill=fill, font=small_font)
                x += d.textlength(exponent, font=small_font)
                i = end + 2
                continue
        next_math = text.find("$", i)
        chunk = text[i:] if next_math == -1 else text[i:next_math]
        d.text((x, y), chunk, fill=fill, font=base_font)
        x += d.textlength(chunk, font=base_font)
        i += len(chunk)
    return x


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


def draw_log2_label(
    d: ImageDraw.ImageDraw,
    xy: Tuple[float, float],
    suffix: str,
    axis_font,
    *,
    scale: int,
) -> float:
    x, y = xy
    sub_font = font(max(1, int(axis_font.size * 0.78)), True)
    prefix_w = d.textlength("log", font=axis_font)
    sub_w = d.textlength("2", font=sub_font)
    d.text((x, y), "log", fill=(25, 25, 25), font=axis_font)
    sub_x = x + prefix_w + 3 * scale
    d.text((sub_x, y + 26 * scale), "2", fill=(25, 25, 25), font=sub_font)
    d.text((sub_x + sub_w + 8 * scale, y), suffix, fill=(25, 25, 25), font=axis_font)
    return prefix_w + sub_w + d.textlength(suffix, font=axis_font) + 11 * scale


def draw_centered_log2_label(
    d: ImageDraw.ImageDraw,
    center_x: float,
    y: float,
    suffix: str,
    axis_font,
    *,
    scale: int,
) -> None:
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    sub_font = font(max(1, int(axis_font.size * 0.78)), True)
    total_w = (
        measure.textlength("log", font=axis_font)
        + measure.textlength("2", font=sub_font)
        + measure.textlength(suffix, font=axis_font)
        + 11 * scale
    )
    draw_log2_label(d, (center_x - total_w / 2, y), suffix, axis_font, scale=scale)


def draw_rotated_log2_label(
    img: Image.Image,
    text_suffix: str,
    center_x: float,
    center_y: float,
    axis_font,
    *,
    scale: int,
) -> None:
    measure = ImageDraw.Draw(Image.new("RGB", (1, 1)))
    sub_font = font(max(1, int(axis_font.size * 0.78)), True)
    total_w = int(
        measure.textlength("log", font=axis_font)
        + measure.textlength("2", font=sub_font)
        + measure.textlength(text_suffix, font=axis_font)
        + 40 * scale
    )
    tmp = Image.new("RGBA", (total_w, 78 * scale), (255, 255, 255, 0))
    td = ImageDraw.Draw(tmp)
    draw_log2_label(td, (0, 0), text_suffix, axis_font, scale=scale)
    bbox = tmp.getbbox()
    if bbox:
        tmp = tmp.crop(bbox)
    rot = tmp.rotate(90, expand=True)
    img.paste(rot, (int(center_x - rot.width / 2), int(center_y - rot.height / 2)), rot)


def y_tick_label(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2g}"


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


def rect_hits_curves(
    rect: Tuple[float, float, float, float],
    polylines: Sequence[Sequence[Tuple[float, float]]],
    pad: float,
) -> bool:
    padded = (rect[0] - pad, rect[1] - pad, rect[2] + pad, rect[3] + pad)
    for pts in polylines:
        for a, b in zip(pts, pts[1:]):
            if segment_hits_rect(a, b, padded):
                return True
    return False


def choose_inset_rect(
    panel_box: Tuple[float, float, float, float],
    polylines: Sequence[Sequence[Tuple[float, float]]],
    *,
    scale: int,
) -> Tuple[float, float, float, float]:
    x0, y0, x1, y1 = panel_box
    panel_w = x1 - x0
    panel_h = y1 - y0
    preferred = (x0 + 0.70 * panel_w, y0 + 0.48 * panel_h)
    best: Optional[Tuple[float, float, float, float]] = None
    best_score = float("inf")
    for factor in [1.0, 0.94, 0.88, 0.82, 0.76, 0.70]:
        inset_w = 0.60 * panel_w * factor
        inset_h = 0.34 * panel_h * factor
        x_candidates = np.linspace(x0 + 0.04 * panel_w, x1 - inset_w - 0.02 * panel_w, 10)
        y_candidates = np.linspace(y0 + 0.08 * panel_h, y1 - inset_h - 0.05 * panel_h, 14)
        for left in x_candidates:
            for top in y_candidates:
                rect = (float(left), float(top), float(left + inset_w), float(top + inset_h))
                if rect_hits_curves(rect, polylines, 9 * scale):
                    continue
                center = ((rect[0] + rect[2]) / 2.0, (rect[1] + rect[3]) / 2.0)
                score = abs(center[0] - preferred[0]) + 0.75 * abs(center[1] - preferred[1])
                if score < best_score:
                    best_score = score
                    best = rect
        if best is not None:
            return best
    return (x1 - 0.42 * panel_w, y0 + 0.30 * panel_h, x1 - 0.02 * panel_w, y0 + 0.56 * panel_h)


def draw_zoom_inset(
    d: ImageDraw.ImageDraw,
    panel: Panel,
    panel_box: Tuple[float, float, float, float],
    polylines: Sequence[Sequence[Tuple[float, float]]],
    tick_font,
    *,
    scale: int,
) -> None:
    if panel.zoom_series is None:
        return
    zoom_series = {
        name: [(x, y) for x, y in panel.zoom_series.get(name, []) if math.isfinite(float(y))]
        for name in ["PAC-alpha(best)", "PRW alpha=16", "PRW alpha=64", "PRW alpha=256"]
    }
    zoom_series = {name: pts for name, pts in zoom_series.items() if pts}
    if not zoom_series:
        return

    ix0, iy0, ix1, iy1 = choose_inset_rect(panel_box, polylines, scale=scale)
    inset_w = ix1 - ix0
    inset_h = iy1 - iy0
    xmin, xmax = 40.0, 80.0
    all_y = [float(y) for pts in zoom_series.values() for _, y in pts]
    ymin, ymax = min(all_y), max(all_y)
    pad = max(0.4, 0.08 * (ymax - ymin))
    ymin = math.floor(ymin - pad)
    ymax = math.ceil(ymax + pad)
    if ymax <= ymin:
        ymax = ymin + 1.0

    def xm(value: float) -> float:
        return ix0 + (value - xmin) / (xmax - xmin) * inset_w

    def ym(value: float) -> float:
        return iy0 + (ymax - value) / (ymax - ymin) * inset_h

    inset_title_font = font(25 * scale, True)
    inset_tick_font = font(24 * scale)
    d.rectangle([ix0, iy0, ix1, iy1], fill=(255, 255, 255), outline=(120, 120, 120), width=2 * scale)
    for tick in [40, 60, 80]:
        px = xm(float(tick))
        d.line([(px, iy0), (px, iy1)], fill=(232, 232, 232), width=1 * scale)
        label = f"{tick:g}"
        tw = d.textlength(label, font=inset_tick_font)
        d.text((px - tw / 2, iy1 - 31 * scale), label, fill=(55, 55, 55), font=inset_tick_font)
    for tick in np.linspace(ymin, ymax, 3):
        py = ym(float(tick))
        d.line([(ix0, py), (ix1, py)], fill=(232, 232, 232), width=1 * scale)
        label = y_tick_label(float(tick))
        d.text((ix0 + 8 * scale, py - 13 * scale), label, fill=(55, 55, 55), font=inset_tick_font)
    d.text((ix0 + 10 * scale, iy0 + 7 * scale), "σ=40-80", fill=(35, 35, 35), font=inset_title_font)

    for name, pts in zoom_series.items():
        color = color_for(panel, name)
        xy = [(xm(x), ym(y)) for x, y in pts]
        for a, b in zip(xy, xy[1:]):
            d.line([a, b], fill=color, width=3 * scale)


def title_lines(title: str) -> List[str]:
    return [title]


def plot_panel(
    img: Image.Image,
    d: ImageDraw.ImageDraw,
    panel: Panel,
    box: Tuple[float, float, float, float],
    title_font,
    tick_font,
    legend_font,
    *,
    scale: int,
) -> None:
    x0, y0, x1, y1 = box
    panel_w = x1 - x0
    panel_h = y1 - y0
    all_x = [x for pts in panel.series.values() for x, _ in pts]
    all_y = [y for pts in panel.series.values() for _, y in pts if y is not None and math.isfinite(float(y))]
    xmin, xmax = min(all_x), max(all_x)
    ymin = 0.0 if panel.y_zero else math.floor(min(all_y) - 1.0)
    ymax = max(all_y)
    if panel.y_zero:
        ymax = 1.12 * (ymax - ymin) + ymin
    else:
        ymax = math.ceil(ymax + 1.0)
    if ymax <= ymin:
        ymax = ymin + 1.0

    def xm(value: float) -> float:
        return x0 + (value - xmin) / (xmax - xmin) * panel_w

    def ym(value: float) -> float:
        return y0 + (ymax - value) / (ymax - ymin) * panel_h

    lines = title_lines(panel.title)
    line_step = 62 * scale
    title_y = y0 - (58 * scale if len(lines) == 1 else 120 * scale)
    for line_idx, line in enumerate(lines):
        draw_rich_text(d, (x0, title_y + line_idx * line_step), line, title_font, scale=scale)
    d.rectangle([x0, y0, x1, y1], fill=(250, 250, 250), outline=(178, 178, 178))

    for tick in panel.x_ticks:
        px = xm(float(tick))
        d.line([(px, y0), (px, y1)], fill=(226, 226, 226))
        label = f"{tick:g}"
        tw = d.textlength(label, font=tick_font)
        d.text((px - tw / 2, y1 + 12 * scale), label, fill=(45, 45, 45), font=tick_font)

    for tick_idx, tick in enumerate(np.linspace(ymin, ymax, 6)):
        py = ym(float(tick))
        d.line([(x0, py), (x1, py)], fill=(226, 226, 226))
        if tick_idx == 0:
            continue
        label = y_tick_label(float(tick))
        tw = d.textlength(label, font=tick_font)
        d.text((x0 - tw - 11 * scale, py - 13 * scale), label, fill=(45, 45, 45), font=tick_font)

    d.line([(x0, y1), (x1, y1)], fill=(30, 30, 30), width=3 * scale)
    d.line([(x0, y0), (x0, y1)], fill=(30, 30, 30), width=3 * scale)

    items = list(panel.series.items())
    draw_order = list(range(len(items)))
    for idx, (name, _) in enumerate(items):
        if name == "Ground truth":
            draw_order = [idx] + [i for i in draw_order if i != idx]
            break
    polylines: List[List[Tuple[float, float]]] = []
    for idx in draw_order:
        name, pts = items[idx]
        color = color_for(panel, name)
        clean = [(float(x), float(y)) for x, y in pts if y is not None and math.isfinite(float(y))]
        xy = [(xm(x), ym(y)) for x, y in clean]
        polylines.append(xy)
        for a, b in zip(xy, xy[1:]):
            d.line([a, b], fill=color, width=5 * scale)
        if panel.draw_markers:
            for px, py in xy:
                d.ellipse(
                    [px - 4 * scale, py - 4 * scale, px + 4 * scale, py + 4 * scale],
                    fill=color,
                    outline="white",
                    width=1 * scale,
                )

    if panel.zoom_series is not None:
        draw_zoom_inset(d, panel, box, polylines, tick_font, scale=scale)

    if panel.inline_legend:
        label_x = x0 + 22 * scale
        label_y = y0 + 20 * scale
        for idx, name in enumerate(panel.series.keys()):
            color = color_for(panel, name)
            yy = label_y + idx * 64 * scale
            d.line(
                [(label_x, yy + 16 * scale), (label_x + 42 * scale, yy + 16 * scale)],
                fill=color,
                width=5 * scale,
            )
            d.text((label_x + 52 * scale, yy), name, fill=(28, 28, 28), font=legend_font)


def draw_legend(
    d: ImageDraw.ImageDraw,
    labels: Sequence[str],
    colors: Dict[str, Tuple[int, int, int]],
    center_x: float,
    y: float,
    legend_font,
    *,
    scale: int,
) -> None:
    widths = [
        52 * scale + d.textlength(display_curve_label(label), font=legend_font) + 28 * scale
        for label in labels
    ]
    x = center_x - sum(widths) / 2.0
    for label, width in zip(labels, widths):
        text = display_curve_label(label)
        color = colors[label]
        x = draw_horizontal_label(d, (x, y), text, color, legend_font, scale=scale)


def render(path: Path, panels: Sequence[Panel]) -> None:
    final_width, final_height = 5600, 2500
    scale = 2
    width, height = final_width * scale, final_height * scale
    img = Image.new("RGB", (width, height), "white")
    d = ImageDraw.Draw(img)

    title_font = font(56 * scale, True)
    tick_font = font(56 * scale)
    axis_font = font(56 * scale, True)
    legend_font = font(56 * scale)

    left = 260 * scale
    right = 60 * scale
    gap_x = 220 * scale
    panel_w = (width - left - right - 3 * gap_x) / 4.0
    panel_h = 850 * scale
    top_y = 190 * scale
    bottom_y = 1450 * scale

    draw_legend(d, GAUSSIAN_ORDER, GAUSSIAN_COLORS, width / 2.0, 38 * scale, legend_font, scale=scale)

    for row, y0 in enumerate([top_y, bottom_y]):
        for col in range(4):
            idx = row * 4 + col
            x0 = left + col * (panel_w + gap_x)
            plot_panel(
                img,
                d,
                panels[idx],
                (x0, y0, x0 + panel_w, y0 + panel_h),
                title_font,
                tick_font,
                legend_font,
                scale=scale,
            )

    top_center = width / 2.0
    top_xlabel = "Gaussian noise std σ"
    tw = d.textlength(top_xlabel, font=axis_font)
    d.text((top_center - tw / 2, 1140 * scale), top_xlabel, fill=(25, 25, 25), font=axis_font)
    draw_rotated_log2_label(
        img,
        " of posterior success rate",
        60 * scale,
        top_y + panel_h / 2,
        axis_font,
        scale=scale,
    )

    f_x0 = left + 1 * (panel_w + gap_x)
    h_x1 = left + 3 * (panel_w + gap_x) + panel_w
    draw_legend(
        d,
        MECHANISM_ORDER,
        MECHANISM_COLORS,
        (f_x0 + h_x1) / 2.0,
        1220 * scale,
        legend_font,
        scale=scale,
    )

    ef_center = left + (2 * panel_w + gap_x) / 2.0
    gh_left = left + 2 * (panel_w + gap_x)
    gh_center = gh_left + (2 * panel_w + gap_x) / 2.0
    draw_centered_log2_label(
        d,
        ef_center,
        2405 * scale,
        " of target posterior success rate",
        axis_font,
        scale=scale,
    )
    iter_label = "Composition iteration"
    tw = d.textlength(iter_label, font=axis_font)
    d.text((gh_center - tw / 2, 2405 * scale), iter_label, fill=(25, 25, 25), font=axis_font)

    for col, label in [(0, "Gaussian noise std σ"), (1, "Expected absolute perturbation")]:
        x0 = left + col * (panel_w + gap_x)
        label_x = x0 - ((200 if col == 0 else 160) * scale)
        draw_rotated_label(
            img,
            label,
            label_x,
            bottom_y + panel_h / 2,
            axis_font,
            scale=scale,
        )

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    img = img.resize((final_width, final_height), resampling)
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def main() -> None:
    panels = gaussian_panels(load_rows(GAUSSIAN_CSV), load_rows(GAUSSIAN_ZOOM_CSV)) + mechanism_panels(load_rows(MECHANISM_CSV))
    render(OUT / "power_experiments_2x4.png", panels)
    render(OUT / "fig1_8_combined_2x4_v11_no_line_overlap_magnifiers.png", panels)
    print(f"wrote {OUT / 'power_experiments_2x4.png'}")
    print(f"wrote {OUT / 'fig1_8_combined_2x4_v11_no_line_overlap_magnifiers.png'}")


if __name__ == "__main__":
    main()
