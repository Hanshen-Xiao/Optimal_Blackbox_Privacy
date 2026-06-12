"""Generate the learning-theoretic privatization workflow figure."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "learning_privacy_workflow.png"


def font(size: int, bold: bool = False):
    candidates = []
    if bold:
        candidates.extend([
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
            "/Library/Fonts/Arial Bold.ttf",
        ])
    candidates.extend([
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
    ])
    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except Exception:
            pass
    return ImageFont.load_default()


def rounded_rect(draw: ImageDraw.ImageDraw, box, radius, fill, outline, width=2):
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def centered_text(draw: ImageDraw.ImageDraw, box, lines, text_font, fill=(25, 25, 25), gap=8):
    x0, y0, x1, y1 = box
    heights = []
    widths = []
    for line in lines:
        bbox = draw.textbbox((0, 0), line, font=text_font)
        widths.append(bbox[2] - bbox[0])
        heights.append(bbox[3] - bbox[1])
    total_h = sum(heights) + gap * (len(lines) - 1)
    y = y0 + (y1 - y0 - total_h) / 2
    for line, w, h in zip(lines, widths, heights):
        draw.text((x0 + (x1 - x0 - w) / 2, y), line, font=text_font, fill=fill)
        y += h + gap


def arrow(draw: ImageDraw.ImageDraw, start, end, color=(70, 70, 70), width=5):
    draw.line([start, end], fill=color, width=width)
    x0, y0 = start
    x1, y1 = end
    size = 16
    if x1 >= x0:
        pts = [(x1, y1), (x1 - size, y1 - size * 0.65), (x1 - size, y1 + size * 0.65)]
    else:
        pts = [(x1, y1), (x1 + size, y1 - size * 0.65), (x1 + size, y1 + size * 0.65)]
    draw.polygon(pts, fill=color)


def main() -> None:
    scale = 2
    w, h = 2400 * scale, 760 * scale
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)

    title_font = font(46 * scale, True)
    box_font = font(30 * scale, True)
    small_font = font(24 * scale)
    label_font = font(25 * scale, True)

    title = "Learning-Theoretic View of Black-Box Privatization"
    tw = d.textlength(title, font=title_font)
    d.text(((w - tw) / 2, 44 * scale), title, font=title_font, fill=(20, 20, 20))

    boxes = [
        ((95, 175, 405, 330), ["Specified prior D", "and black-box leakage F"], (235, 246, 255)),
        ((500, 175, 810, 330), ["Sample secrets", "x_i ~ D"], (242, 249, 236)),
        ((905, 175, 1215, 330), ["Evaluate leakage", "F(x_i, v)"], (242, 249, 236)),
        ((1310, 175, 1675, 330), ["Empirical PRW", "optimization"], (255, 245, 230)),
        ((1770, 175, 2205, 330), ["Minimum-utility", "empirical mechanism"], (255, 245, 230)),
    ]

    lower_boxes = [
        ((515, 485, 875, 635), ["Validation /", "confidence correction"], (248, 240, 255)),
        ((1015, 485, 1425, 635), ["Population posterior-risk", "guarantee"], (235, 246, 255)),
        ((1580, 485, 2115, 635), ["Asymptotically tight", "population optimum"], (235, 246, 255)),
    ]

    for box, lines, fill in boxes + lower_boxes:
        box = tuple(int(v * scale) for v in box)
        rounded_rect(d, box, 18 * scale, fill, (120, 130, 140), 3 * scale)
        centered_text(d, box, lines, box_font)

    for i in range(len(boxes) - 1):
        b0 = boxes[i][0]
        b1 = boxes[i + 1][0]
        arrow(
            d,
            (int(b0[2] * scale), int((b0[1] + b0[3]) / 2 * scale)),
            (int(b1[0] * scale), int((b1[1] + b1[3]) / 2 * scale)),
        )

    arrow(d, (int(1985 * scale), int(330 * scale)), (int(695 * scale), int(485 * scale)))
    arrow(d, (int(875 * scale), int(560 * scale)), (int(1015 * scale), int(560 * scale)))
    arrow(d, (int(1425 * scale), int(560 * scale)), (int(1580 * scale), int(560 * scale)))

    notes = [
        ((1320, 365), "Optimization gap: eta", (180, 95, 20)),
        ((890, 672), "Statistical gap: validation/generalization", (95, 55, 150)),
        ((1545, 672), "Finite-alpha gap vanishes as alpha -> infinity", (30, 110, 90)),
    ]
    for (x, y), text, color in notes:
        d.text((x * scale, y * scale), text, font=label_font, fill=color)

    caption = "No white-box structure of F is required: the privacy design problem is learned from sampled secrets and black-box evaluations."
    wrapped = [
        "No white-box structure of F is required: the privacy design problem is learned from",
        "sampled secrets and black-box evaluations."
    ]
    y = 695 * scale
    for line in wrapped:
        lw = d.textlength(line, font=small_font)
        d.text(((w - lw) / 2, y), line, font=small_font, fill=(65, 65, 65))
        y += 32 * scale

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    img = img.resize((w // scale, h // scale), resampling)
    img.save(OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
