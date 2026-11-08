#!/usr/bin/env python3
"""Create a readable Word copy of the source-backed dataset guide."""
from pathlib import Path
import subprocess

from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

DATA = Path(__file__).resolve().parents[1]
SOURCES = DATA / 'dataset_guide_sources_20261106'
SOURCE = DATA / '数据集说明_来源与测量内容.md'
TARGET = SOURCE.with_suffix('.docx')

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21), Cm(29.7)
sec.top_margin = sec.bottom_margin = sec.left_margin = sec.right_margin = Cm(2)
for name in ['Normal', 'Body Text', 'First Paragraph', 'Compact']:
    try:
        s = doc.styles[name]
    except KeyError:
        s = doc.styles.add_style(name, 1)
    s.font.name = 'Noto Sans CJK SC'
    s.font.size = Pt(10.5)
    s._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), 'Noto Sans CJK SC')
    s.paragraph_format.line_spacing = 1.25
    s.paragraph_format.space_after = Pt(6)
for name, size in [('Title', 21), ('Heading 1', 18), ('Heading 2', 14), ('Heading 3', 11.5)]:
    s = doc.styles[name]
    s.font.name = 'Noto Sans CJK SC'
    s.font.size = Pt(size)
    s.font.bold = True
    s.font.color.rgb = RGBColor.from_string('17365D')
    s._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), 'Noto Sans CJK SC')
    s.paragraph_format.space_before = Pt(12)
    s.paragraph_format.space_after = Pt(6)
    s.paragraph_format.keep_with_next = True
footer = sec.footer.paragraphs[0]
footer.alignment = 1
field = OxmlElement('w:fldSimple')
field.set(qn('w:instr'), 'PAGE')
footer._p.append(field)
reference = SOURCES / 'reference.docx'
doc.save(reference)

subprocess.run(['pandoc', str(SOURCE), '--from=gfm', '--to=docx', '--standalone',
                '--reference-doc=' + str(reference), '--resource-path=' + str(DATA),
                '--metadata=lang:zh-CN', '--output=' + str(TARGET)], check=True)
doc = Document(TARGET)
for table in doc.tables:
    header = OxmlElement('w:tblHeader')
    table.rows[0]._tr.get_or_add_trPr().append(header)
    for row in table.rows:
        keep = OxmlElement('w:cantSplit')
        row._tr.get_or_add_trPr().append(keep)
        for cell in row.cells:
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(3)
                for run in p.runs:
                    run.font.size = Pt(8.5 if len(row.cells) >= 4 else 10)
                    run.font.name = 'Noto Sans CJK SC'
doc.save(TARGET)
print(TARGET)
