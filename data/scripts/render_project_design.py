#!/usr/bin/env python3
"""Render the established local project design to a readable Word companion."""
from pathlib import Path
import re
import subprocess
from docx import Document
from docx.shared import Cm, Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '项目总设计_20261106'


def main():
    names = ['README.md', 'data_audit.md', 'integration_design.md', 'model_design.md', 'experiment_plan.md', 'literature_and_reporting.md']
    parts = []
    for name in names:
        text = (OUT / name).read_text()
        text = re.sub(r'```mermaid\n.*?```', '流程：原始数据与协议 → 身份/条件/对照核查 → 主实验、异构监督及外部表达 → 基线/共享模型 → 分组泛化与早期信息检验 → 按证据组织主张。', text, flags=re.S)
        parts.append(text)
    combined = OUT / '合订本来源.md'
    combined.write_text('\n\n'.join(parts))
    target = OUT / '项目整体设计与数据对接_20261106.docx'
    reference = ROOT / 'data/dataset_guide_sources_20261106/reference.docx'
    subprocess.run(['pandoc', str(combined), '--from=gfm', '--to=docx', '--standalone', '--toc',
                    '--metadata=lang:zh-CN', '--reference-doc=' + str(reference), '--output=' + str(target)], check=True)
    d = Document(target)
    for s in d.sections:
        s.left_margin = s.right_margin = Cm(1.8)
    for table in d.tables:
        h = OxmlElement('w:tblHeader')
        table.rows[0]._tr.get_or_add_trPr().append(h)
        for row in table.rows:
            keep = OxmlElement('w:cantSplit')
            row._tr.get_or_add_trPr().append(keep)
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    paragraph.paragraph_format.space_after = Pt(3)
                    for run in paragraph.runs:
                        run.font.name = 'Noto Sans CJK SC'
                        run.font.size = Pt(8 if len(row.cells) > 3 else 9)
                        run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'), 'Noto Sans CJK SC')
    d.save(target)
    print(f'{target}: {len(d.paragraphs)} paragraphs, {len(d.tables)} tables')


if __name__ == '__main__':
    main()
