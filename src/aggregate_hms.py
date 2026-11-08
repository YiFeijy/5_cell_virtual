#!/usr/bin/env python3
"""Recompute HMS ratios from the three original files and check preserved labels."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def condition_id(row):
    value = f'{row.library_id}|{row.cell_id}|{row.pert_id}|{row.dose_uM:.9g}'
    return 'HMSv2_LIB_' + hashlib.sha256(value.encode()).hexdigest()[:16]


def main():
    frames = []
    for path in sorted((ROOT/'data/raw/HMS_LINCS').glob('*.xls')):
        d = pd.read_excel(path, sheet_name='Data')
        d['source_file'] = path.name
        d['original_excel_row'] = np.arange(len(d))+2
        frames.append(d)
    assert len(frames) == 3
    raw = pd.concat(frames, ignore_index=True).rename(columns={
        'CellLine':'cell_id','BRDid':'pert_id','Small Mol Concentration':'dose_uM',
        'Relative Cell Count':'relative_cell_count','Relative Growth Rate After Treatment':'relative_net_growth'})
    assert raw.Timepoint.eq(72).all() and raw['Small Mol Conc Unit'].str.lower().eq('um').all()
    raw['library_id'] = raw['LJP Library Plate'].astype(str).str.extract(r'(\d+)')[0].astype(int).map(lambda x:f'LJP{x:03d}')
    raw['dose_uM'] = pd.to_numeric(raw.dose_uM)
    keys = ['library_id','cell_id','pert_id','dose_uM']
    bio = raw.groupby(keys+['source_file']).agg(relative_cell_count=('relative_cell_count','mean'),
                                               relative_net_growth=('relative_net_growth','mean')).reset_index()
    means = bio.groupby(keys).relative_cell_count.mean().reset_index()
    means['library_condition_id'] = [condition_id(r) for r in means.itertuples()]
    original = pd.read_csv(ROOT/'data/reference/hms_conditions.csv')
    check = means.merge(original[['library_condition_id','relative_cell_count']], on='library_condition_id', suffixes=('_recomputed','_preserved'), validate='one_to_one')
    assert len(check) == len(original) == 4140
    error = float(abs(check.relative_cell_count_recomputed-check.relative_cell_count_preserved).max())
    assert error < 1e-12
    print(json.dumps({'original_wells':len(raw),'biological_repeat_records':len(bio),
                      'library_conditions':len(check),'max_ratio_difference':error,
                      'labels_checked':'passed','files_modified':False}, indent=2))


if __name__ == '__main__':
    main()
