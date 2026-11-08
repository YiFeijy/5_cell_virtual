#!/usr/bin/env python3
"""Trace a large cross-library discrepancy back to actual HMS wells."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / '数据可行性调研_20261105/evidence/xlrd-2.0.2-py2.py3-none-any.whl'))
OUT = ROOT / 'data/investigation_20261106/CGP60474_跨库板差异'


def main():
    OUT.mkdir(exist_ok=True)
    parts = []
    for p in sorted((ROOT / 'data/raw/HMS_LINCS').glob('*.xls')):
        d = pd.read_excel(p, sheet_name='Data')
        d['source_excel_row'] = np.arange(len(d)) + 2
        d['source_file'] = p.name
        d = d[d.BRDid.eq('BRD-K79090631')].copy()
        d['recomputed_count_ratio'] = d['Total Cell Count After Treatment'] / d['Total Control Cell Count']
        parts.append(d)
    w = pd.concat(parts, ignore_index=True)
    # Biofile 2 records the two synonyms "CGP-60474;CGP60474".
    assert w.DrugName.str.upper().str.split(';').apply(lambda names: all(n.replace('-', '') == 'CGP60474' for n in names)).all()
    assert w.Timepoint.eq(72).all() and w['Small Mol Conc Unit'].str.lower().eq('um').all()
    difference = abs(w['recomputed_count_ratio'] - w['Relative Cell Count'])
    assert difference.max() < .000051
    w.to_csv(OUT / 'CGP60474_全部原孔及重算.csv', index=False)
    b = w.groupby(['CellLine', 'LJP Library Plate', 'Small Mol Concentration', 'source_file']).agg(
        count_ratio=('Relative Cell Count', 'mean'), n_wells=('source_excel_row', 'size')).reset_index()
    b.to_csv(OUT / '库板内生物重复.csv', index=False)
    c = b.groupby(['CellLine', 'LJP Library Plate', 'Small Mol Concentration']).agg(
        count_ratio=('count_ratio', 'mean'), sd=('count_ratio', 'std'), n_bio=('source_file', 'size')).reset_index()
    c.to_csv(OUT / '两库板各自剂量曲线.csv', index=False)
    case = b[b.CellLine.eq('SKBR3') & b['Small Mol Concentration'].eq(.04)]
    result = {'drug': 'CGP-60474', 'pert_id': 'BRD-K79090631', 'raw_CGP_wells': len(w),
              'nominal_dose_and_units_confirmed': True, 'count_formula_max_absolute_error': float(difference.max()),
              'SKBR3_0.04uM_bioreplicates': case.to_dict('records'),
              'interpretation': 'The library difference is present in each source biological file and survives direct Nt/Nc recomputation. It is not created by cross-library averaging or this project\'s ratio calculation. Cause (library-specific material, actual exposure, protocol or source annotation) remains unresolved.'}
    (OUT / 'audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
