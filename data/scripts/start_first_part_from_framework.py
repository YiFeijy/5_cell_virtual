#!/usr/bin/env python3
"""Read the user's framework without editing it and audit local state references."""
from pathlib import Path
import hashlib
import json
import re
import zipfile
from collections import Counter
import h5py
import pandas as pd
from docx import Document
from docx.text.paragraph import Paragraph

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '第一部分_依据整体框架的设计与起步_20261106'
SOURCE = ROOT / '面向抗癌药物筛选的虚拟细胞建模研究.docx'


def main():
    OUT.mkdir(exist_ok=True)
    d = Document(SOURCE)
    records = []
    part = ''
    for i, element in enumerate(d.element.body):
        if element.tag.split('}')[-1] != 'p':
            continue
        p = Paragraph(element, d)
        if p.style.name == 'Heading 2':
            part = p.text
        if p.text.strip():
            records.append({'body_block_index': i, 'style': p.style.name, 'research_part': part, 'text': p.text})
    frame = pd.DataFrame(records)
    frame.to_csv(OUT / '整体框架原文提取.csv', index=False)
    first = frame[frame.research_part.str.startswith('1、')]
    first.to_csv(OUT / '第一部分原文提取.csv', index=False)
    assert len(first) == 6 and len(d.tables) == 0
    source = {'path': str(SOURCE.relative_to(ROOT)), 'sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'bytes': SOURCE.stat().st_size, 'paragraph_count_including_blank': len(d.paragraphs),
              'table_count': len(d.tables), 'inline_images': len(d.inline_shapes),
              'first_part_body_blocks': first.body_block_index.tolist(), 'source_modified_by_this_script': False}
    (OUT / '框架来源核对.json').write_text(json.dumps(source, ensure_ascii=False, indent=2) + '\n')

    observations = []
    pattern = re.compile(r'viab|live|dead|death|annexin', re.I)
    for path in sorted((ROOT / 'data/raw/sciPlex').glob('*.h5ad')):
        with h5py.File(path) as f:
            fields = list(f['obs'])
            observations.append({'dataset': path.name, 'path': str(path.relative_to(ROOT)),
                                 'shape': list(map(int, f['X'].attrs['shape'])), 'obs_fields': fields,
                                 'explicit_live_dead_like_field_names': [k for k in fields if pattern.search(k)],
                                 'interpretation': 'Field-name audit only. Captured expression and QC annotations do not by themselves authenticate independently measured viable-cell state. Original protocol and supplementary references require matching.'})
    counts = Counter()
    schemas = Counter()
    zip_count = 0
    for path in sorted((ROOT / 'data/raw/MIXseq').glob('*.zip')):
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                if name.startswith('__MACOSX/') or not name.endswith('/classifications.csv'):
                    continue
                meta = pd.read_csv(archive.open(name))
                counts.update(meta.cell_quality.fillna('MISSING').astype(str))
                schemas.update(meta.columns)
                zip_count += 1
    observations.append({'dataset': 'MIXseq_downloaded_archives', 'classification_tables_checked': zip_count,
                         'cell_quality_value_counts': dict(counts), 'metadata_field_occurrence_counts': dict(schemas),
                         'explicit_live_dead_like_field_names': [k for k in schemas if pattern.search(k)],
                         'interpretation': 'Author cell_quality labels retained. The label normal is not relabeled as a measured live-cell gate; quality definition/protocol must be traced separately.'})
    (OUT / '单细胞状态参考字段核查.json').write_text(json.dumps(observations, ensure_ascii=False, indent=2) + '\n')
    rows = [
        ['HMS计数', 'N0/Nt/Ncontrol；72h数量与净增长', '独立于表达的细胞数量读数', '相对数量不是存活比例，也不是亚群组成', '主功能监督及重复可靠性'],
        ['LINCS L5', '3/24h群体扰动响应z', '响应目标；非独立活细胞状态参考', '不凭计数比例倒推残存细胞表达', '响应解码与共享表征监督'],
        ['LINCS L3与vehicle', '同平台处理/同步对照表达', '表达尺度一致的参考', '同步对照非同一处理样本真实t0；参考QC待做', '状态重建与drug-control差异敏感性方案'],
        ['MIXseq', '捕获的单细胞表达、细胞身份和QC标签', '可研究捕获群体的组成与组内变化', 'normal标签非已认证活细胞门；附带sens源于GDSC/PRISM', '组成分层/伪bulk/时间与表达外部评价'],
        ['sciPlex2/3', '单细胞/核表达及条件注释', '可研究捕获群体分布', '本地obs字段未发现明确live/dead字段；独立存活参考需追协议/附件', '组成与表达分布评价；有对应独立参考后开展存活相关估计']
    ]
    pd.DataFrame(rows, columns=['source', 'actual_observation', 'reference_support', 'question_to_resolve', 'use_now']).to_csv(OUT / '状态观测与参考支持表.csv', index=False)
    mapping = [
        [13, '研究总目标', '数据关联、表达与功能关系、联合预测和药后状态表征', '全部第一部分方案'],
        [14, '多来源关联', '同项目主队列＋可比条件关联＋未配对任务监督', '数据层及条件误差比较'],
        [15, '状态及组成', '群体响应起步；单细胞分层；有独立参考时估计存活相关状态', '状态定义与参考子研究'],
        [16, '联合模型', '背景/结构/剂量/时间→受两任务监督的状态表示→两个观测输出', '最小联合模型及部分标签扩展'],
        [17, '独立评价与筛选', '单任务/常规多任务对照、未见药/条件、状态和排序评价', '比较实验与应用评价']
    ]
    pd.DataFrame(mapping, columns=['source_body_block_index', 'framework_requirement', 'implementation', 'design_location']).to_csv(OUT / '第一部分逐项对应.csv', index=False)
    print(json.dumps({'source': source, 'state_reference_sources_checked': len(observations), 'mixseq_classification_tables': zip_count}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
