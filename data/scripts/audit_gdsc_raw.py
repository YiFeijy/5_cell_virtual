#!/usr/bin/env python3
"""Count actual GDSC well records, assays, endpoint durations and controls."""
from pathlib import Path
from collections import Counter
import json
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1]
OUT = DATA / 'analysis_20261106'
results = {}
for family in ['GDSC1', 'GDSC2']:
    p = next((DATA / 'unpacked' / family).rglob('*.csv'))
    n = treated = bad_treatment = 0
    duration, assays, tags = Counter(), Counter(), Counter()
    cell_ids, drug_ids, plates, nc_plates, blank_plates = set(), set(), set(), set(), set()
    for d in pd.read_csv(p, usecols=['BARCODE', 'CELL_LINE_NAME', 'DRUG_ID', 'CONC', 'INTENSITY', 'TAG', 'DURATION', 'ASSAY'], chunksize=250000):
        n += len(d)
        duration.update(d.DURATION.fillna('MISSING').value_counts().to_dict())
        assays.update(d.ASSAY.fillna('MISSING').value_counts().to_dict())
        tags.update(d.TAG.fillna('MISSING').value_counts().to_dict())
        cell_ids.update(d.CELL_LINE_NAME.dropna().unique())
        plates.update(d.BARCODE.unique())
        nc_plates.update(d.loc[d.TAG.isin(['NC-0', 'NC-1']), 'BARCODE'].unique())
        blank_plates.update(d.loc[d.TAG.eq('B'), 'BARCODE'].unique())
        t = d[d.TAG.fillna('').str.match(r'^L\d+-D\d+-S$')]
        treated += len(t)
        drug_ids.update(t.DRUG_ID.dropna().unique())
        bad_treatment += int((~np.isfinite(t.INTENSITY) | ~np.isfinite(t.CONC) | t.CONC.le(0)).sum())
    results[family] = {'raw_rows': n, 'cell_names': len(cell_ids), 'plates': len(plates),
        'library_single_treatment_rows': treated, 'library_drug_ids': len(drug_ids),
        'library_rows_with_invalid_intensity_or_dose': bad_treatment,
        'duration_days_raw_row_counts': dict(duration), 'assay_code_raw_row_counts': dict(assays),
        'plates_with_NC0_or_NC1': len(nc_plates), 'plates_with_blank': len(blank_plates),
        'plates_without_negative_control': len(plates - nc_plates),
        'qualification': 'Raw intensity and tags, not normalized viability. DURATION is days. Raw rows include controls and potentially repeated tags; not a count of independent treatments.'}
    pd.DataFrame(tags.most_common(), columns=['TAG', 'rows']).to_csv(OUT / (family + '_well_tag_counts.csv'), index=False)
    print(family, json.dumps(results[family]), flush=True)
    (OUT / 'GDSC_actual_raw_audit.json').write_text(json.dumps(results, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))
