"""Build the review PDF from rGEDI's generated Rd documentation."""
from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, KeepTogether, PageBreak, PageTemplate,
    Paragraph, Preformatted, Spacer, Table, TableStyle
)
from reportlab.platypus.tableofcontents import TableOfContents

ROOT = Path(__file__).resolve().parents[1]
TEXT = ROOT / "output" / ".reference-text"
OUTPUT = ROOT / "output" / "pdf" / "rGEDI-reference-manual.pdf"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

meta = {}
for line in (TEXT / "metadata.tsv").read_text(encoding="utf-8").splitlines():
    key, value = line.split("\t", 1)
    meta[key] = value

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="CoverTitle", parent=styles["Title"], fontName="Helvetica-Bold",
                          fontSize=30, leading=34, textColor=colors.HexColor("#123B2A"),
                          alignment=TA_CENTER, spaceAfter=15))
styles.add(ParagraphStyle(name="CoverSub", parent=styles["Normal"], fontSize=13, leading=18,
                          textColor=colors.HexColor("#315B4A"), alignment=TA_CENTER))
styles.add(ParagraphStyle(name="Topic", parent=styles["Heading1"], fontSize=18, leading=22,
                          textColor=colors.HexColor("#123B2A"), spaceAfter=10))
styles.add(ParagraphStyle(name="Section", parent=styles["Heading2"], fontSize=10, leading=13,
                          textColor=colors.HexColor("#28754E"), spaceBefore=7, spaceAfter=3))
styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], fontSize=8.3,
                          leading=11, spaceAfter=4))
code_style = ParagraphStyle(name="Code", fontName="Courier", fontSize=6.9, leading=8.6,
                            leftIndent=8, rightIndent=4, textColor=colors.HexColor("#17352A"),
                            backColor=colors.HexColor("#F2F6F3"), borderPadding=5)

class ManualDoc(BaseDocTemplate):
    def __init__(self, filename):
        super().__init__(filename, pagesize=letter, rightMargin=0.62*inch,
                         leftMargin=0.62*inch, topMargin=0.62*inch,
                         bottomMargin=0.58*inch, title="rGEDI Function Reference")
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="normal")
        self.addPageTemplates(PageTemplate(id="main", frames=frame, onPage=self._page))

    def _page(self, canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#D5E4DB"))
        canvas.line(self.leftMargin, 0.43*inch, letter[0]-self.rightMargin, 0.43*inch)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(colors.HexColor("#567066"))
        canvas.drawString(self.leftMargin, 0.25*inch, f"rGEDI {meta.get('Version', '')} function reference")
        canvas.drawRightString(letter[0]-self.rightMargin, 0.25*inch, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph) and flowable.style.name == "Topic":
            title = flowable.getPlainText()
            key = "topic-%s" % re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(title, key, level=0, closed=True)
            self.notify("TOCEntry", (0, title, self.page, key))

def esc(text):
    return (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))

def split_sections(text):
    text = text.replace("\b", "")
    lines = text.splitlines()
    sections, current, bucket = [], "OVERVIEW", []
    known = {"NAME", "DESCRIPTION", "USAGE", "ARGUMENTS", "DETAILS", "VALUE",
             "REFERENCES", "AUTHOR(S)", "SEE ALSO", "EXAMPLES", "NOTE", "FORMAT"}
    for line in lines:
        candidate = line.strip()
        heading = candidate.rstrip(":").upper()
        if heading in known:
            if any(x.strip() for x in bucket):
                sections.append((current, "\n".join(bucket).strip()))
            current, bucket = heading, []
        else:
            bucket.append(line.rstrip())
    if any(x.strip() for x in bucket):
        sections.append((current, "\n".join(bucket).strip()))
    return sections

story = [Spacer(1, 0.38*inch), Paragraph("rGEDI", styles["CoverTitle"]),
         Paragraph("GEDI data access, processing, simulation, modeling, and mapping",
                   styles["CoverSub"]), Spacer(1, 0.28*inch)]
cover = ROOT / "readme" / "fig-gedi-wall-to-wall.png"
if cover.exists():
    story += [Image(str(cover), width=6.45*inch, height=3.64*inch), Spacer(1, 0.24*inch)]
story += [Paragraph(f"Function reference • version {meta.get('Version', '')} • {meta.get('Date', '')}",
                    styles["CoverSub"]),
          Paragraph(f"{meta.get('Documents', '')} documented help topics generated from the package source",
                    styles["CoverSub"]), PageBreak(),
          Paragraph("Contents", styles["Topic"])]
toc = TableOfContents()
toc.levelStyles = [ParagraphStyle(name="TOC0", fontName="Helvetica", fontSize=8.3,
                                  leading=11, leftIndent=8, firstLineIndent=-8,
                                  textColor=colors.HexColor("#244B3B"))]
story += [toc, PageBreak(), Paragraph("Workflow at a glance", styles["Topic"]),
          Paragraph("The README examples proceed from NASA CMR search and authenticated download "
                    "through local or cloud reading, extraction, quality filtering, clipping, "
                    "rasterization, modeling, simulation, and Earth Engine wall-to-wall mapping.",
                    styles["BodySmall"]), Spacer(1, 8)]
gallery = []
for name in ["fig-study-site.png", "fig-waveform-rh.png", "fig-simulator.png", "fig-gedi-model-validation.png"]:
    path = ROOT / "readme" / name
    if path.exists():
        gallery.append(Image(str(path), width=3.1*inch, height=2.12*inch))
if len(gallery) == 4:
    table = Table([[gallery[0], gallery[1]], [gallery[2], gallery[3]]], colWidths=[3.25*inch]*2)
    table.setStyle(TableStyle([("VALIGN", (0,0), (-1,-1), "MIDDLE"),
                               ("BOX", (0,0), (-1,-1), 0.4, colors.HexColor("#C9DDD1")),
                               ("INNERGRID", (0,0), (-1,-1), 0.4, colors.HexColor("#C9DDD1")),
                               ("BACKGROUND", (0,0), (-1,-1), colors.white)]))
    story += [table]
story += [PageBreak()]

for path in sorted(TEXT.glob("*.txt"), key=lambda p: p.stem.lower()):
    raw = path.read_text(encoding="utf-8", errors="replace")
    sections = split_sections(raw)
    title = path.stem
    story.append(Paragraph(esc(title), styles["Topic"]))
    for heading, body in sections:
        if heading in {"OVERVIEW", "NAME"} or not body:
            continue
        story.append(Paragraph(heading.title(), styles["Section"]))
        if heading in {"USAGE", "EXAMPLES"}:
            story.append(Preformatted(body, code_style, maxLineLength=105))
        else:
            paragraphs = [p.strip() for p in re.split(r"\n\s*\n", body) if p.strip()]
            for paragraph in paragraphs:
                story.append(Paragraph(esc(re.sub(r"\s+", " ", paragraph)), styles["BodySmall"]))
    story.append(PageBreak())

doc = ManualDoc(str(OUTPUT))
doc.multiBuild(story)
print(OUTPUT)
