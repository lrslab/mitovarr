"""Draw the README overview as editable, journal-sized vector artwork.

Run with Python and ReportLab; add --png for a 600-dpi Poppler rendering.
The five panels are schematics, not observations or benchmark measurements.
This documentation utility is independent of the Rust build and runtime.
"""

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Drawing, Group, Line, Path as VPath
from reportlab.graphics.shapes import Polygon, Rect, String
from reportlab.lib.colors import HexColor, white
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas


OUT = Path(__file__).resolve().parent
WIDTH, HEIGHT = 183, 168
INK = HexColor("#202124")
GRAY = HexColor("#6D7680")
RULE = HexColor("#C5CBD0")
LIGHT = HexColor("#E7E9EB")
PAPER = HexColor("#F5F6F7")
BLUE = HexColor("#0072B2")
PALE_BLUE = HexColor("#E1EFF7")
ORANGE = HexColor("#D55E00")
PALE_ORANGE = HexColor("#F9E8DC")
GREEN = HexColor("#009E73")
PALE_GREEN = HexColor("#E0F1EA")


def register_fonts():
    directory = Path(os.environ.get(
        "MITOVARR_FIGURE_FONT_DIR", "/System/Library/Fonts/Supplemental"
    ))
    for name, filename in (("Arial", "Arial.ttf"), ("Arial-Bold", "Arial Bold.ttf")):
        source = directory / filename
        if not source.is_file():
            raise SystemExit(
                f"Missing {source}. Set MITOVARR_FIGURE_FONT_DIR to the directory "
                "containing Arial.ttf and Arial Bold.ttf."
            )
        pdfmetrics.registerFont(TTFont(name, str(source)))


def text(group, x, y, value, size=6.3, bold=False, anchor="start"):
    """Coordinates are millimetres from the upper left; y is the baseline."""
    font = "Arial-Bold" if bold else "Arial"
    width = pdfmetrics.stringWidth(value, font, size) / mm
    left = x - {"start": 0, "middle": width / 2, "end": width}[anchor]
    if not (0 <= left and left + width <= WIDTH and 0 < y < HEIGHT):
        raise ValueError(f"Text outside artwork: {value!r}")
    group.add(String(x * mm, (HEIGHT - y) * mm, value, fontName=font,
                     fontSize=size, textAnchor=anchor, fillColor=INK))


def line(group, x1, y1, x2, y2, color=GRAY, width=.55, dash=None):
    group.add(Line(x1 * mm, (HEIGHT - y1) * mm, x2 * mm, (HEIGHT - y2) * mm,
                   strokeColor=color, strokeWidth=width, strokeDashArray=dash))


def polyline(group, points, color=GRAY, width=.65, dash=None):
    path = VPath(fillColor=None, strokeColor=color, strokeWidth=width,
                 strokeDashArray=dash)
    for index, (x, y) in enumerate(points):
        if index == 0:
            path.moveTo(x * mm, (HEIGHT - y) * mm)
        else:
            path.lineTo(x * mm, (HEIGHT - y) * mm)
    group.add(path)


def rect(group, x, y, w, h, fill=None, stroke=None, width=.5, dash=None):
    group.add(Rect(x * mm, (HEIGHT - y - h) * mm, w * mm, h * mm,
                   fillColor=fill, strokeColor=stroke, strokeWidth=width,
                   strokeDashArray=dash))


def arrow(group, x1, y1, x2, y2, color=GRAY, width=.6, head=1.15, dash=None):
    line(group, x1, y1, x2, y2, color, width, dash)
    angle = math.atan2(y2 - y1, x2 - x1)
    points = [(x2, y2),
              (x2 - head * math.cos(angle) + .42 * head * math.sin(angle),
               y2 - head * math.sin(angle) - .42 * head * math.cos(angle)),
              (x2 - head * math.cos(angle) - .42 * head * math.sin(angle),
               y2 - head * math.sin(angle) + .42 * head * math.cos(angle))]
    group.add(Polygon([v for x, y in points for v in (x * mm, (HEIGHT - y) * mm)],
                      fillColor=color, strokeColor=None))


def pair(group, left, right, y, read=4, color=BLUE):
    # The fine line joins mates from one fragment; it is not sequenced coverage.
    line(group, left, y, right, y, GRAY, .35)
    arrow(group, left, y, left + read, y, color, 1, .8)
    arrow(group, right, y, right - read, y, color, 1, .8)


def block(group, x, y, w, label=None, fill=PALE_BLUE, stroke=BLUE,
          height=3.8, dashed=False):
    rect(group, x, y, w, height, fill, stroke, .6,
         [1.5, 1.2] if dashed else None)
    if label:
        text(group, x + w / 2, y + height / 2 + .72, label, 6, anchor="middle")


def arc(group, cx, cy, radius, start, end, color=GRAY, width=.7):
    path = VPath(fillColor=None, strokeColor=color, strokeWidth=width)
    count = max(1, math.ceil(abs(end - start) / 45))
    for index in range(count):
        a = math.radians(start + (end - start) * index / count)
        b = math.radians(start + (end - start) * (index + 1) / count)
        f = 4 / 3 * math.tan((b - a) / 4)
        x0, y0 = cx + radius * math.cos(a), cy + radius * math.sin(a)
        x3, y3 = cx + radius * math.cos(b), cy + radius * math.sin(b)
        if index == 0:
            path.moveTo(x0 * mm, (HEIGHT - y0) * mm)
        path.curveTo((x0 - f * radius * math.sin(a)) * mm,
                     (HEIGHT - y0 - f * radius * math.cos(a)) * mm,
                     (x3 + f * radius * math.sin(b)) * mm,
                     (HEIGHT - y3 + f * radius * math.cos(b)) * mm,
                     x3 * mm, (HEIGHT - y3) * mm)
    group.add(path)


def panel(drawing, letter, x, y, title):
    group = Group()
    drawing.add(group, name=f"panel_{letter}")
    text(group, x, y, letter, 8, True)
    text(group, x + 4.8, y, title, 7, True)
    return group


def gene_arc(group, cx, cy, radius, start, end):
    """A filled gene segment; outlines keep publication-scale line weights."""
    points = []
    for r, angles in ((radius + .45, range(start, end + 1)),
                      (radius - .45, range(end, start - 1, -1))):
        points.extend((cx + r * math.cos(math.radians(a)),
                       cy + r * math.sin(math.radians(a))) for a in angles)
    group.add(Polygon([v for x, y in points for v in (x * mm, (HEIGHT - y) * mm)],
                      fillColor=GREEN, strokeColor=None))


def make_drawing():
    drawing = Drawing(WIDTH * mm, HEIGHT * mm,
                      initialFontName="Arial", initialFontSize=6.3)
    drawing.add(Rect(0, 0, WIDTH * mm, HEIGHT * mm,
                     fillColor=white, strokeColor=None))

    # a: input routes retain their distinct assembly and validation scopes.
    a = panel(drawing, "a", 3, 5.5, "Input routes define the available evidence")
    for left, right, color in ((5, 59, BLUE), (64, 118, BLUE), (127, 179, GRAY)):
        line(a, left, 8.5, right, 8.5, color, .8)
    text(a, 32, 12.5, "BAM-baited", 6.7, True, "middle")
    text(a, 32, 16.8, "Indexed mitochondrial BAM", 6.1, anchor="middle")
    pair(a, 19, 38, 21.2, read=3.6)
    pair(a, 28, 45, 23.8, read=3.6)
    text(a, 32, 28, "MT query + mate recovery", 6, anchor="middle")
    text(a, 32, 32, "Validate extracted pairs", 6, anchor="middle")

    text(a, 91, 12.5, "Seeded", 6.7, True, "middle")
    text(a, 91, 16.8, "Original pairs + seed FASTA", 6.1, anchor="middle")
    block(a, 75, 20.2, 9, "Seed", fill=LIGHT, stroke=GRAY, height=3.4)
    arrow(a, 86, 21.9, 92, 21.9)
    pair(a, 94, 107, 21.9, read=3)
    text(a, 91, 28, "Recruit complete pairs", 6, anchor="middle")
    text(a, 91, 32, "Validate all original pairs", 6, anchor="middle")

    text(a, 153, 12.5, "Mapping assembly", 6.7, True, "middle")
    text(a, 153, 16.8, "Reference + indexed BAM", 6.1, anchor="middle")
    line(a, 139, 23, 171, 23, GRAY, .7)
    pair(a, 140, 156, 20.5, read=3)
    pair(a, 154, 170, 22, read=3)
    text(a, 153, 28, "Fragment-aware pileup", 6, anchor="middle")
    text(a, 153, 32, "Local repair: matching FASTQ", 6, anchor="middle")
    for center in (32, 91):
        line(a, center, 33.5, center, 35, BLUE)
    line(a, 32, 35, 91, 35, BLUE)
    arrow(a, 45, 35, 45, 37, BLUE, head=.8)
    arrow(a, 153, 33.5, 153, 37, GRAY, head=.8)

    # b: repeat topology is proposed from sample reads, not reference length.
    b = panel(drawing, "b", 3, 41, "Preserve alternative structures in repeat-rich regions")
    text(b, 5, 47, "De novo: sample-read graph", 6.5, True)
    text(b, 100, 47, "Mapping: local sample-read reconstruction", 6.5, True)
    line(b, 92, 46, 92, 71.5, RULE, .4)

    for x, label in ((7, "A"), (23, "R"), (39, "B")):
        repeat = label == "R"
        block(b, x, 58, 6.5, label, height=3.8,
              fill=PALE_ORANGE if repeat else PALE_BLUE,
              stroke=ORANGE if repeat else BLUE)
    arrow(b, 13.5, 59.9, 23, 59.9, BLUE)
    arrow(b, 29.5, 59.9, 39, 59.9, BLUE)
    polyline(b, [(26.3, 58), (26.3, 52.5), (35, 52.5),
                 (35, 56), (30, 56)], ORANGE)
    arrow(b, 30, 56, 27.8, 58, ORANGE, head=.85)
    arrow(b, 46.8, 59.9, 52.5, 59.9)
    line(b, 54, 55, 54, 65)
    line(b, 52.5, 59.9, 54, 59.9)
    for y, labels in ((53.1, ["A", "R", "B"]),
                       (63.1, ["A", "R", "R", "B"])):
        arrow(b, 54, y + 1.9, 56, y + 1.9, head=.7)
        line(b, 56, y + 1.9, 56 + 7.7 * len(labels) - 1.2, y + 1.9)
        for index, label in enumerate(labels):
            repeat = label == "R"
            block(b, 56 + index * 7.7, y, 6.5, label,
                  fill=PALE_ORANGE if repeat else PALE_BLUE,
                  stroke=ORANGE if repeat else BLUE)
    text(b, 5, 71.2, "Read bridges + bounded shorter-k connections", 6.1)

    text(b, 100, 55.4, "Pileup", 6)
    text(b, 100, 65.4, "Proposal", 6)
    for y in (53, 63):
        line(b, 115, y + 1.9, 177, y + 1.9)
        block(b, 115, y, 7, "A", fill=LIGHT, stroke=GRAY)
        block(b, 130, y, 7, "R", fill=PALE_ORANGE, stroke=ORANGE)
        block(b, 170, y, 7, "B", fill=LIGHT, stroke=GRAY)
    rect(b, 143, 52.6, 21, 4.6, white)
    line(b, 143, 54.9, 164, 54.9, GRAY, .55, [1.5, 1.4])
    text(b, 153.5, 52.3, "uncertain", 6, anchor="middle")
    block(b, 143, 63, 21, "Sample path", fill=PALE_ORANGE, stroke=ORANGE)
    for x in (141, 166):
        line(b, x, 50.5, x, 68.5, RULE, .4, [1.4, 1.4])
    text(b, 100, 71.2, "Whole-contig tests retain homologous source copies", 6.1)
    for x in (45, 148):
        line(b, x, 72.8, x, 74.5)
    line(b, 45, 74.5, 148, 74.5)
    arrow(b, 45, 74.5, 45, 77, head=.85)

    # c: fragment likelihood and copy-depth fit are different evidence uses.
    c = panel(drawing, "c", 3, 81, "Read evidence and depth inference")
    text(c, 5, 86, "Compare the same fragments across candidates", 6.3, True)
    for left, right, y in ((7, 30, 90), (12, 34, 94), (8, 32, 98)):
        pair(c, left, right, y, read=4.2)
    text(c, 40, 90.7, "Gapped read alignments", 6.2)
    text(c, 40, 94.7, "Mate geometry, if estimable", 6.2)
    text(c, 40, 98.7, "Observed junction support", 6.2)
    line(c, 5, 102, 87, 102, RULE, .45)
    text(c, 5, 106.5, "Optional multi-k / coverage layouts", 6.3, True)
    # Explicitly illustrative normalized depth, not a data-derived chart.
    text(c, 20, 110.6, "Relative R1 k-mer depth", 6, anchor="middle")
    line(c, 8, 112, 8, 122, GRAY, .4)
    line(c, 8, 122, 32, 122, GRAY, .4)
    for value, y in ((0, 122), (1, 118), (2, 114)):
        line(c, 7.3, y, 8, y, GRAY, .4)
        text(c, 6.2, y + .7, str(value), 6, anchor="end")
    for x, y, label, fill, stroke in ((11, 118, "A", PALE_BLUE, BLUE),
                                     (19, 114, "R", PALE_ORANGE, ORANGE),
                                     (27, 118, "B", PALE_BLUE, BLUE)):
        rect(c, x, y, 4, 122 - y, fill, stroke, .45)
        text(c, x + 2, 125.5, label, 6, anchor="middle")
    text(c, 20, 129.4, "Illustrative", 5.8, anchor="middle")
    text(c, 39, 113.5, "Pool sample paths across k", 6.1)
    text(c, 39, 117.5, "Fit copy usage to unique flanks", 6.1)
    text(c, 39, 121.5, "Retain competing read evidence", 6.1)
    text(c, 39, 127, "Depth does not establish order", 6.1, True)
    text(c, 5, 134, "Search limits and unscored evidence remain explicit", 6)

    # d: three expected output states, with no simulated performance claims.
    d = panel(drawing, "d", 99, 81, "Expected outputs")
    for y, h, color, fill, dash in ((85, 17.5, BLUE, PALE_BLUE, None),
                                   (105, 12, GRAY, PAPER, None),
                                   (119.5, 13, ORANGE, PALE_ORANGE, [2, 1.4])):
        rect(d, 100, y, 79, h, fill, color, .55, dash)
    text(d, 103, 89, "Strict checks passed", 6.5, True)
    text(d, 103, 93.3, "Final consensus + evidence", 6.2)
    text(d, 103, 97.3, "De novo: read-supported circle", 6)
    text(d, 103, 100.7, "Mapping: circularity unverified", 6)
    for x in (160, 166, 172):
        block(d, x, 89.7, 4.5, height=2.7, fill=PALE_BLUE, stroke=BLUE)
    line(d, 158, 91, 160, 91, BLUE)
    line(d, 164.5, 91, 166, 91, BLUE)
    line(d, 170.5, 91, 172, 91, BLUE)
    text(d, 103, 109, "Partial or unresolved candidates", 6.5, True)
    text(d, 103, 112.7, "Partial sequence + uncertainty report", 6)
    text(d, 103, 115.9, "Retain alternatives; withhold final", 6)
    polyline(d, [(162, 110), (166, 110), (169, 107.8), (176, 107.8)], GRAY, .6)
    polyline(d, [(166, 110), (169, 112.2), (176, 112.2)], GRAY, .6)
    text(d, 103, 123.4, "Depth-inferred layout", 6.5, True)
    text(d, 103, 127.2, "Sequence hypothesis + copy-depth fit", 6)
    text(d, 103, 130.6, "Unverified; strict final status unchanged", 6)
    for x in (162, 169):
        block(d, x, 122, 5, height=2.8, fill=PALE_ORANGE,
              stroke=ORANGE, dashed=True)
    line(d, 167, 123.4, 169, 123.4, ORANGE, .5, [1, 1])
    # A strict comparison may resolve or remain unresolved. The optional
    # depth-derived branch reaches only the unverified-layout output.
    line(d, 88.5, 94, 94, 94)
    arrow(d, 94, 94, 98.5, 94, head=1)
    line(d, 94, 94, 94, 111)
    arrow(d, 94, 111, 98.5, 111, head=1)
    arrow(d, 88.5, 125, 98.5, 125, ORANGE, .55, 1, [2, 1.4])
    text(d, 100, 137.1, "Freeze sequence + status", 6)
    arrow(d, 144, 134.4, 144, 138.7, head=1)

    # e: annotation is independent of graph inference and assembly selection.
    e = panel(drawing, "e", 3, 144, "Biological interpretation after sequence selection")
    text(e, 179, 144, "Species gene models enter here", 6, anchor="end")
    line(e, 5, 153, 56, 153)
    block(e, 5, 151.1, 18, "CDS")
    block(e, 26, 151.1, 12, "tRNA", fill=PALE_GREEN, stroke=GREEN)
    block(e, 41, 151.1, 15, "rRNA", fill=LIGHT, stroke=GRAY)
    text(e, 5, 161.2, "Homology + gene-quality checks", 6.1)
    for cx, start, end in ((75, 5, 40), (106, 270, 305)):
        arc(e, cx, 153.5, 4.3, 0, 360, width=.65)
        gene_arc(e, cx, 153.5, 4.3, start, end)
    arrow(e, 82.5, 153.5, 98.2, 153.5)
    text(e, 112, 151, "1", 6)
    text(e, 92, 161.2, "tRNA-Pro origin, when eligible", 6.1, anchor="middle")
    line(e, 129, 153, 179, 153)
    block(e, 129, 151.1, 10, "Gene", fill=LIGHT, stroke=GRAY)
    block(e, 143, 151.1, 21, "NCR", fill=white, stroke=GRAY)
    block(e, 169, 151.1, 10, "Gene", fill=LIGHT, stroke=GRAY)
    polyline(e, [(145, 150), (145, 148.7), (162, 148.7), (162, 150)],
             ORANGE, .55, [1.5, 1.3])
    text(e, 129, 161.2, "Positional NCRs + uncertainty", 6.1)

    key = Group()
    drawing.add(key, name="figure_key")
    text(key, 5, 166.3,
         "Schematic examples, not measured results. A/B, unique flanks; R, repeat. Annotation preserves assembly status.",
         5.8)
    return drawing


def export(drawing, png, dpi):
    stem = OUT / "package_workflow"
    canvas = Canvas(str(stem.with_suffix(".pdf")),
                    pagesize=(WIDTH * mm, HEIGHT * mm),
                    initialFontName="Arial", initialFontSize=6.3,
                    pageCompression=1, invariant=1)
    canvas.setTitle("Figure 1 | Evidence-aware reconstruction of repeat-rich mitochondrial genomes with mitovarr")
    canvas.setAuthor("mitovarr")
    canvas.setSubject("Five schematic panels; editable RGB vector artwork, 183 x 168 mm")
    renderPDF.draw(drawing, canvas, 0, 0)
    canvas.showPage()
    canvas.save()
    renderSVG.drawToFile(drawing, str(stem.with_suffix(".svg")))
    namespace = "http://www.w3.org/2000/svg"
    ET.register_namespace("", namespace)
    tree = ET.parse(stem.with_suffix(".svg"))
    root = tree.getroot()
    root.set("width", f"{WIDTH}mm")
    root.set("height", f"{HEIGHT}mm")
    root.set("viewBox", f"0 0 {WIDTH * mm:.6f} {HEIGHT * mm:.6f}")
    root.set("role", "img")
    root.set("aria-labelledby", "title description")
    for tag, value in (("title", "Figure 1. Evidence-aware reconstruction of repeat-rich mitochondrial genomes with mitovarr"),
                       ("desc", "Five schematic panels: distinct input and validation scopes; sample-derived alternative repeat structures; shared-fragment evidence and optional multi-k copy-depth inference; conditional consensus, partial or unresolved candidates and unverified layout hypotheses; independent gene annotation, eligible tRNA-Pro normalization and positional NCRs. All sequences and the relative-depth example are illustrative, not experimental results.")):
        element = ET.Element(f"{{{namespace}}}{tag}", id="description" if tag == "desc" else tag)
        element.text = value
        root.insert(0, element)
    wrapper = root.find(f"{{{namespace}}}g")
    wrapper.set("transform", f"scale(1,-1) translate(0,{-HEIGHT * mm:.6f})")
    content = wrapper.find(f"{{{namespace}}}g")
    groups = [node for node in content if node.tag == f"{{{namespace}}}g"]
    names = ["panel_a", "panel_b", "panel_c", "panel_d", "panel_e", "figure_key"]
    if len(groups) != len(names):
        raise ValueError("Unexpected SVG groups; cannot label editable panels")
    for name, group in zip(names, groups):
        group.set("id", name)
    # ReportLab's internal bold font name is not a portable CSS family.
    # Also retain the exact baseline rather than its rounded SVG reflection.
    for node in root.findall(f".//{{{namespace}}}text"):
        style = node.get("style", "")
        if "font-family: Arial-Bold;" in style:
            style = style.replace("font-family: Arial-Bold;",
                                  "font-family: Arial, Helvetica, sans-serif; font-weight: 700;")
        else:
            style = style.replace("font-family: Arial;",
                                  "font-family: Arial, Helvetica, sans-serif; font-weight: 400;")
        node.set("style", style)
        node.set("transform", f"translate(0,{2 * float(node.get('y')):.8f}) scale(1,-1)")
    tree.write(stem.with_suffix(".svg"), encoding="utf-8", xml_declaration=True)
    if png:
        renderer = shutil.which("pdftoppm")
        if not renderer:
            raise SystemExit("PDF/SVG saved. Put pdftoppm on PATH to export PNG.")
        subprocess.run([renderer, "-singlefile", "-r", str(dpi), "-png",
                        str(stem.with_suffix(".pdf")), str(stem)], check=True)
    print(f"Saved PDF/SVG at {WIDTH} x {HEIGHT} mm" +
          (f" and {dpi}-dpi PNG" if png else ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png", action="store_true", help="Render PNG with pdftoppm")
    parser.add_argument("--dpi", type=int, default=600, help="PNG dpi (default: 600)")
    args = parser.parse_args()
    if not 72 <= args.dpi <= 1200:
        parser.error("--dpi must be between 72 and 1200")
    register_fonts()
    export(make_drawing(), args.png, args.dpi)


if __name__ == "__main__":
    main()
