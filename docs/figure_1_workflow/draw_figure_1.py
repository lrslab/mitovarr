"""Draw the current mitovarr workflow as editable publication artwork.

Usage: python draw_figure_1.py --png --dpi 600
Requires ReportLab; PNG export additionally uses Poppler's pdftoppm.
This documentation utility is independent of the Rust build and runtime.
All depicted reads, graphs and loci are schematic, not experimental results.
"""

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Circle, Drawing, Group, Line, Path as VPath, Polygon, Rect, String
from reportlab.lib.colors import HexColor, white
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas


OUT = Path(__file__).resolve().parent
WIDTH, HEIGHT = 183, 205
# Compact vertical spacing without scaling type or changing the column width.
Y_SCALE = HEIGHT / 237
INK = HexColor("#202124")
GRAY = HexColor("#7B838A")
LIGHT = HexColor("#DDE1E5")
RULE = HexColor("#C3C9CE")
PAPER = HexColor("#F5F6F7")
BLUE = HexColor("#0072B2")
PALE_BLUE = HexColor("#EAF3F8")
ORANGE = HexColor("#D55E00")
PALE_ORANGE = HexColor("#FBF0E9")
PURPLE = HexColor("#CC79A7")
PALE_PURPLE = HexColor("#F4E5EE")
GREEN = HexColor("#009E73")

fonts = Path(os.environ.get("MITOVARR_FIGURE_FONT_DIR", "/System/Library/Fonts/Supplemental"))
for name, file in (("Arial", "Arial.ttf"), ("Arial-Bold", "Arial Bold.ttf")):
    pdfmetrics.registerFont(TTFont(name, str(fonts / file)))

drawing = Drawing(WIDTH * mm, HEIGHT * mm, initialFontName="Arial", initialFontSize=6.3)
drawing.add(Rect(0, 0, WIDTH * mm, HEIGHT * mm, fillColor=white, strokeColor=None))
panels = []


def text(g, x, y, value, size=6.3, bold=False, anchor="start"):
    """Use design coordinates from the upper left; y is the text baseline."""
    g.add(String(x * mm, (HEIGHT - y * Y_SCALE) * mm, value,
                 fontName="Arial-Bold" if bold else "Arial", fontSize=size,
                 fillColor=INK, textAnchor=anchor))


def line(g, x1, y1, x2, y2, color=GRAY, width=.55, dash=None):
    g.add(Line(x1 * mm, (HEIGHT - y1 * Y_SCALE) * mm, x2 * mm, (HEIGHT - y2 * Y_SCALE) * mm,
               strokeColor=color, strokeWidth=width, strokeDashArray=dash))


def rect(g, x, y, w, h, fill=None, stroke=None, width=.5, dash=None):
    g.add(Rect(x * mm, (HEIGHT - (y + h) * Y_SCALE) * mm, w * mm, h * Y_SCALE * mm,
               fillColor=fill, strokeColor=stroke, strokeWidth=width, strokeDashArray=dash))


def arrow(g, x1, y1, x2, y2, color=GRAY, width=.65, head=1.15):
    line(g, x1, y1, x2, y2, color, width)
    angle = math.atan2(y2 - y1, x2 - x1)
    points = [(x2, y2),
              (x2 - head * math.cos(angle) + .42 * head * math.sin(angle),
               y2 - head * math.sin(angle) - .42 * head * math.cos(angle)),
              (x2 - head * math.cos(angle) - .42 * head * math.sin(angle),
               y2 - head * math.sin(angle) + .42 * head * math.cos(angle))]
    g.add(Polygon([v for x, y in points for v in (x * mm, (HEIGHT - y * Y_SCALE) * mm)],
                  fillColor=color, strokeColor=None))


def pair(g, x1, x2, y, color=BLUE, read=4, size=.9):
    """Inward read pair; the connecting line denotes one DNA fragment."""
    line(g, x1, y, x2, y, GRAY, .35)
    arrow(g, x1, y, x1 + read, y, color, size, .9)
    arrow(g, x2, y, x2 - read, y, color, size, .9)


def panel(letter, x, y, title):
    g = Group()
    drawing.add(g)
    panels.append((letter, g))
    text(g, x, y, letter, 8, True)
    text(g, x + 5, y, title, 7, True)
    return g


def block(g, x, y, w, h=3.4, fill=LIGHT, label=None, stroke=None):
    rect(g, x, y, w, h, fill, stroke, .45)
    if label:
        text(g, x + w / 2, y + h / 2 + .75, label, 6, anchor="middle")


def arc(g, cx, cy, radius, start, end, color=GRAY, width=.8):
    # Keep circular molecules circular, despite the tighter vertical rhythm.
    cy *= Y_SCALE
    path = VPath(fillColor=None, strokeColor=color, strokeWidth=width)
    n = max(1, int(math.ceil(abs(end - start) / 45)))
    for i in range(n):
        a = math.radians(start + (end - start) * i / n)
        b = math.radians(start + (end - start) * (i + 1) / n)
        f = 4 / 3 * math.tan((b - a) / 4)
        x0, y0 = cx + radius * math.cos(a), cy + radius * math.sin(a)
        x3, y3 = cx + radius * math.cos(b), cy + radius * math.sin(b)
        if i == 0:
            path.moveTo(x0 * mm, (HEIGHT - y0) * mm)
        path.curveTo((x0 - f * radius * math.sin(a)) * mm,
                     (HEIGHT - y0 - f * radius * math.cos(a)) * mm,
                     (x3 + f * radius * math.sin(b)) * mm,
                     (HEIGHT - y3 + f * radius * math.cos(b)) * mm,
                     x3 * mm, (HEIGHT - y3) * mm)
    g.add(path)


# a: BAM-baited and seeded de novo inputs, alongside the BAM mapping route.
a = panel("a", 4, 5.7, "Input routes and read-evidence scope")
for x in (5, 64):
    rect(a, x, 10, 54, 33, PAPER)
rect(a, 127, 10, 51, 33, PAPER)
text(a, 32, 15.2, "BAM-baited", 6.7, True, "middle")
text(a, 32, 20, "CaenDR mtDNA BAM", 6.2, anchor="middle")
text(a, 32, 24.1, "MT query + mate recovery", 5.9, anchor="middle")
arrow(a, 32, 25.4, 32, 29.1, BLUE)
text(a, 32, 32.5, "Extracted read pairs", 6.3, True, "middle")
text(a, 32, 39.2, "Validate: extracted pairs", 5.8, anchor="middle")

text(a, 91, 15.2, "Seeded", 6.7, True, "middle")
text(a, 91, 20, "Original pairs + mt seed", 6.2, anchor="middle")
text(a, 91, 24.1, "Iterative pair recruitment", 5.9, anchor="middle")
arrow(a, 91, 25.4, 91, 29.1, BLUE)
text(a, 91, 32.5, "Recruited read pairs", 6.3, True, "middle")
text(a, 91, 39.2, "Validate: all original pairs", 5.8, anchor="middle")

text(a, 152.5, 15.2, "Mapping assembly", 6.7, True, "middle")
text(a, 152.5, 20, "Matching reference FASTA", 6.2, anchor="middle")
text(a, 152.5, 24.1, "+ sorted, indexed BAM", 6.2, anchor="middle")
line(a, 135, 34, 170, 34, GRAY, 1)
for x, y in ((137, 28), (146, 30.3), (154, 32.5)):
    pair(a, x, x + 14, y, ORANGE, read=4.1, size=.8)
text(a, 152.5, 39.2, "CaenDR mt.bam or local BAM", 5.8, anchor="middle")

for x in (32, 91):
    line(a, x, 43, x, 45.8, BLUE, .6)
line(a, 32, 45.8, 91, 45.8, BLUE, .6)
text(a, 61.5, 49.7, "Sample FASTQ only; no reference in the assembler", 6.1, anchor="middle")
arrow(a, 46, 50.8, 46, 53.6, BLUE)
arrow(a, 152.5, 43, 152.5, 53.6, ORANGE)

# b: the native graph engine shared by the two active FASTQ entry modes.
b = panel("b", 4, 57, "Sample-graph de novo assembly")
line(b, 5, 60.4, 87, 60.4, BLUE, 1.4)
text(b, 5, 67, "Build, simplify and compact the sample graph", 6.5, True)

nodes = {"a": (9, 79), "b": (16, 79), "c": (23, 74), "d": (23, 84),
         "e": (31, 74), "f": (31, 84), "g": (38, 79), "h": (45, 79)}
for u, v in (("a", "b"), ("b", "c"), ("b", "d"), ("c", "e"),
             ("d", "f"), ("e", "g"), ("f", "g"), ("g", "h")):
    line(b, *nodes[u], *nodes[v], BLUE, .8)
line(b, 16, 79, 16, 71, LIGHT, .8)
for x, y in nodes.values():
    b.add(Circle(x * mm, (HEIGHT - y * Y_SCALE) * mm, .85 * mm,
                 fillColor=white, strokeColor=BLUE, strokeWidth=.75))
b.add(Circle(16 * mm, (HEIGHT - 71 * Y_SCALE) * mm, .8 * mm,
             fillColor=white, strokeColor=LIGHT, strokeWidth=.6))
arrow(b, 48, 79, 55, 79, GRAY)
for yy in (74.5, 83.5):
    line(b, 63, 79, 68, yy, BLUE, .7)
    line(b, 74, yy, 79, 79, BLUE, .7)
block(b, 57, 77.3, 6, fill=PALE_BLUE, stroke=BLUE)
block(b, 68, 72.8, 6, fill=PALE_PURPLE, label="R", stroke=PURPLE)
block(b, 68, 81.8, 6, fill=PALE_BLUE, stroke=BLUE)
block(b, 79, 77.3, 6, fill=PALE_BLUE, stroke=BLUE)
text(b, 26.5, 90, "Local error-path pruning", 6.1, anchor="middle")
text(b, 71, 90, "Unitigs", 6.1, anchor="middle")
arrow(b, 46, 92, 46, 97, BLUE)

text(b, 5, 102, "Thread reads and propose complete paths", 6.5, True)
for y, name, two in ((109.2, "P1", False), (121.4, "P2", True)):
    text(b, 8, y + 2.5, name, 6.2)
    line(b, 18, y + 1.7, 84, y + 1.7, BLUE, .6)
    block(b, 18, y, 15, fill=PALE_BLUE, stroke=BLUE)
    block(b, 69, y, 15, fill=PALE_BLUE, stroke=BLUE)
    if two:
        block(b, 39, y, 10, fill=PALE_PURPLE, label="R", stroke=PURPLE)
        block(b, 52, y, 10, fill=PALE_PURPLE, label="R", stroke=PURPLE)
    else:
        block(b, 45.5, y, 10, fill=PALE_PURPLE, label="R", stroke=PURPLE)
pair(b, 26, 76, 116.2, BLUE, read=9)
text(b, 5, 132.6, "Read threading + smaller-k end connections", 6.1)
text(b, 5, 136.5, "Bounded repeat-copy and circular-path hypotheses", 6.1)
rect(b, 5, 140, 82, 17, PALE_BLUE)
text(b, 46, 145.1, "Compare paths using paired fragments", 6.5, True, "middle")
text(b, 46, 150, "Gapped alignments + empirical mate geometry", 6.1, anchor="middle")
text(b, 46, 154.1, "Shared read cohort; each fragment contributes once", 5.9, anchor="middle")
text(b, 46, 162, "Primary sequence + alternatives", 6.3, True, "middle")
text(b, 46, 166, "Final: resolved path, coverage and supported joins", 5.9, anchor="middle")

# c: reference-coordinate consensus with an optional sample-graph extension.
c = panel("c", 95, 57, "Reference-guided mapping assembly")
line(c, 96, 60.4, 178, 60.4, ORANGE, 1.4)
text(c, 96, 67, "Call a fragment-aware pileup consensus", 6.5, True)
for x1, x2, y in ((102, 136, 72), (114, 153, 75.2), (129, 170, 78.4)):
    pair(c, x1, x2, y, ORANGE, read=8.3)
rect(c, 143, 73.9, 1.1, 2.6, ORANGE)
rect(c, 143, 77.1, 1.1, 2.6, ORANGE)
line(c, 101, 84.2, 173, 84.2, GRAY, 1.6)
rect(c, 142.5, 82.6, 2, 3.2, PALE_ORANGE, ORANGE, .5)
text(c, 137, 90, "Quality-weighted SNP / indel calls; low support = N", 6, anchor="middle")
arrow(c, 137, 92, 137, 97, ORANGE)
line(c, 137, 93.5, 181, 93.5, GRAY, .45)
line(c, 181, 93.5, 181, 159, GRAY, .45)
arrow(c, 181, 159, 137, 159, GRAY, .45, .8)
text(c, 163, 96.5, "Pileup only", 5.7, anchor="middle")

rect(c, 96, 98, 82, 60, None, RULE, .6, [2, 2])
text(c, 137, 102, "Optional local assembly with matching FASTQ", 6.3, True, "middle")
text(c, 137, 107.6, "Detect clips, indels and ambiguous regions", 6.1, anchor="middle")
for y, label in ((113.1, "Backbone"), (123.4, "Sample path")):
    text(c, 99, y + 2.4, label, 5.8)
    line(c, 119, y + 1.6, 174, y + 1.6, GRAY, .55)
    block(c, 119, y, 10, fill=LIGHT)
    block(c, 164, y, 10, fill=LIGHT)
block(c, 141.7, 112.8, 9.6, h=4, fill=white)
line(c, 142, 114.8, 151, 114.8, GRAY, .6, [1.5, 1.5])
block(c, 135, 123.4, 23, fill=PALE_PURPLE, label="NCR", stroke=PURPLE)
line(c, 129, 111.5, 129, 129, RULE, .4, [1.5, 1.5])
line(c, 164, 111.5, 164, 129, RULE, .4, [1.5, 1.5])
text(c, 137, 132.6, "Local sample paths between unique flanks", 6.1, anchor="middle")
text(c, 137, 136.5, "Retain homologous copies in whole-contig tests", 6.1, anchor="middle")
rect(c, 98, 140, 78, 17, PALE_ORANGE)
text(c, 137, 145.1, "Compare paths and both junctions", 6.5, True, "middle")
text(c, 137, 150, "Same paired fragments for all candidates", 6.1, anchor="middle")
text(c, 137, 154.1, "Patch supported sequence; preserve uncertainty", 5.9, anchor="middle")
text(c, 137, 162, "Backbone + accepted patches / unresolved candidates", 6.1, True, "middle")
text(c, 137, 166, "Mapping completion does not verify circularity", 5.9, anchor="middle")

# d: annotation is a downstream command; models cannot choose an assembly.
d = panel("d", 4, 175, "Post-assembly annotation")
text(d, 178, 175, "Species homology models enter only here", 6.1, anchor="end")
for x in (46, 137):
    line(d, x, 167.5, x, 170, GRAY, .55)
line(d, 46, 170, 137, 170, GRAY, .55)
arrow(d, 91.5, 170, 91.5, 179, GRAY, .65)
rect(d, 5, 179, 173, 10.3, PAPER)
text(d, 91.5, 183.1, "Frozen FASTA + selection status + assembly ambiguity", 6.5, True, "middle")
text(d, 91.5, 187.2, "No annotation feedback to assembly", 6, anchor="middle")
arrow(d, 31, 189.3, 31, 193, GRAY, .5, .9)

text(d, 5, 197, "Homology-based gene annotation", 6.4, True)
text(d, 66, 197, "Optional tRNA-Pro rotation", 6.4, True)
text(d, 127, 197, "Positional NCR1-4 annotation", 6.4, True)

line(d, 7, 204.5, 55, 204.5, GRAY, .5)
block(d, 7, 202.6, 17, h=3.8, fill=PALE_BLUE, label="CDS", stroke=BLUE)
block(d, 27, 202.6, 9, h=3.8, fill=HexColor("#E2F3ED"), label="tRNA", stroke=GREEN)
block(d, 39, 202.6, 16, h=3.8, fill=PALE_ORANGE, label="rRNA", stroke=ORANGE)
text(d, 5, 213, "12 CDS / 22 tRNA / 2 rRNA models", 6.1)
text(d, 5, 217, "Boundary + coding checks; retain partials", 5.9)

for cx, start, end in ((79, -15, 18), (108, -90, -57)):
    arc(d, cx, 205, 4.5, 0, 360, GRAY, .85)
    arc(d, cx, 205, 4.5, start, end, GREEN, 2.6)
arrow(d, 89, 205, 98, 205, GRAY, .65)
text(d, 111, 200, "1", 5.8)
text(d, 91.5, 213, "5' tRNA-Pro becomes position 1", 6.1, anchor="middle")
text(d, 91.5, 217, "Circular input; unique, complete Pro", 5.9, anchor="middle")

line(d, 128, 204.5, 177, 204.5, GRAY, .5)
block(d, 128, 202.6, 10, h=3.8, fill=LIGHT, label="Gene")
block(d, 141, 202.6, 22, h=3.8, fill=PALE_PURPLE, label="NCR", stroke=PURPLE)
block(d, 167, 202.6, 10, h=3.8, fill=LIGHT, label="Gene")
text(d, 127, 213, "Intervals from observed gene flanks", 6.1)
text(d, 127, 217, "Sequence + boundaries + path evidence", 5.9)
arrow(d, 58, 204.5, 63, 204.5, GRAY, .55, .95)
arrow(d, 120, 204.5, 125, 204.5, GRAY, .55, .95)

line(d, 5, 222, 178, 222, RULE, .5)
text(d, 5, 227, "Reusable outputs", 6.4, True)
text(d, 40, 227, "Sequences (FASTA)  |  Graphs (GFA)  |  Features (GFF3)  |  Evidence (TSV / JSON)", 6.2)
for x, col, label in ((5, BLUE, "Sample graph"), (42, ORANGE, "Reference-guided"),
                      (86, PURPLE, "Repeat / NCR"), (125, GRAY, "Backbone / data flow")):
    rect(d, x, 231.6, 3, 1.5, col)
    text(d, x + 4.5, 233.1, label, 5.7)


def export():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png", action="store_true", help="Also render a PNG using pdftoppm")
    parser.add_argument("--dpi", type=int, default=600, help="PNG resolution (default: 600)")
    args = parser.parse_args()
    stem = OUT / "figure_1_workflow"
    canvas = Canvas(str(stem.with_suffix(".pdf")), pagesize=(WIDTH * mm, HEIGHT * mm),
                    initialFontName="Arial", initialFontSize=6.3, pageCompression=1)
    canvas.setTitle("Detailed method schematic | Native mitochondrial assembly and NCR interpretation")
    canvas.setAuthor("mitovarr")
    canvas.setCreator("mitovarr documentation / ReportLab")
    canvas.setSubject("Schematic of the implemented workflow, 2026-09-22")
    renderPDF.draw(drawing, canvas, 0, 0)
    canvas.showPage()
    canvas.save()
    renderSVG.drawToFile(drawing, str(stem.with_suffix(".svg")))

    # Preserve editable text and expose panel groups in vector editors.
    svg = "http://www.w3.org/2000/svg"
    inkscape = "http://www.inkscape.org/namespaces/inkscape"
    ET.register_namespace("", svg)
    ET.register_namespace("inkscape", inkscape)
    tree = ET.parse(stem.with_suffix(".svg"))
    root = tree.getroot()
    root.set("width", f"{WIDTH}mm")
    root.set("height", f"{HEIGHT}mm")
    root.set("viewBox", f"0 0 {WIDTH * mm:.6f} {HEIGHT * mm:.6f}")
    title = ET.Element(f"{{{svg}}}title")
    title.text = "Detailed method schematic. Native mitochondrial assembly and NCR interpretation with mitovarr"
    root.insert(0, title)
    desc = ET.Element(f"{{{svg}}}desc")
    desc.text = ("Four schematic panels: input contracts; native de novo assembly; "
                 "reference-guided local assembly; independent gene annotation, tRNA-Pro rotation "
                 "and NCR analysis. No experimental measurements are depicted.")
    root.insert(1, desc)
    wrapper = root.find(f"{{{svg}}}g")
    wrapper.set("transform", f"scale(1,-1) translate(0,{-HEIGHT * mm:.6f})")
    content = wrapper.find(f"{{{svg}}}g")
    groups = [n for n in content if n.tag == f"{{{svg}}}g"]
    if len(groups) != len(panels):
        raise ValueError("Unexpected SVG structure: cannot identify the four editable panels")
    for (letter, _), group in zip(panels, groups):
        group.set("id", f"panel_{letter}")
        group.set(f"{{{inkscape}}}groupmode", "layer")
        group.set(f"{{{inkscape}}}label", f"{letter} - workflow")
    tree.write(stem.with_suffix(".svg"), encoding="utf-8", xml_declaration=True)
    if args.png:
        poppler = shutil.which("pdftoppm")
        if not poppler:
            raise SystemExit("PDF/SVG saved. Install Poppler or put pdftoppm on PATH for PNG export.")
        subprocess.run([poppler, "-singlefile", "-r", str(args.dpi), "-png",
                        str(stem.with_suffix(".pdf")), str(stem)], check=True)
    print(f"Created editable PDF/SVG at {WIDTH} x {HEIGHT} mm" +
          (f" and {args.dpi}-dpi PNG" if args.png else ""))


if __name__ == "__main__":
    export()
