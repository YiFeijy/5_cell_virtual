#!/usr/bin/env python3
"""Inspect GEO sample annotations without loading the full microarray matrix."""
from pathlib import Path
import csv
import gzip
import json
import pandas as pd

DATA = Path(__file__).resolve().parents[1]
OUT = DATA / 'analysis_20261106'
rows = {}
with gzip.open(DATA / 'raw/NCI_TPW/GSE116436_series_matrix.txt.gz', 'rt') as f:
    for line in f:
        a = next(csv.reader([line], delimiter='\t'))
        if not a: continue
        key = a[0]
        if key in ['!Sample_title', '!Sample_geo_accession', '!Sample_source_name_ch1', '!Sample_platform_id']:
            rows[key] = a[1:]
        if key == '!series_matrix_table_begin':
            header = next(csv.reader([next(f)], delimiter='\t'))
            first = next(csv.reader([next(f)], delimiter='\t'))
            break
m = pd.DataFrame({k.removeprefix('!Sample_'): v for k, v in rows.items()})
doses = m.title.str.extract(r'_([\d.]+)nM_([\d.]+)h$')
m['dose_nM'] = pd.to_numeric(doses[0], errors='coerce')
m['time_h'] = pd.to_numeric(doses[1], errors='coerce')
m.to_csv(OUT / 'NCI_sample_metadata.csv', index=False)
r = {'matrix_sample_columns': len(header) - 1, 'first_probe_id': first[0],
    'platforms': m.platform_id.unique().tolist(), 'distinct_source_labels': m.source_name_ch1.nunique(),
    'time_h_counts': m.time_h.value_counts(dropna=False).to_dict(),
    'zero_dose_arrays': int(m.dose_nM.eq(0).sum()), 'title_unparsed_dose_or_time': int(m.time_h.isna().sum()),
    'qualification': 'GEO processed microarray matrix; probe-to-gene mapping not completed; zero-dose controls may be named per compound and need identity/batch audit.'}
(OUT / 'NCI_actual_header_audit.json').write_text(json.dumps(r, indent=2))
print(json.dumps(r))
