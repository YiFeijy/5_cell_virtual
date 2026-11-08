#!/usr/bin/env python3
"""Check cell metadata and Matrix Market dimensions in every downloaded MIX-seq ZIP."""
from pathlib import Path
import json
import zipfile
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1]
OUT = DATA / 'analysis_20261106'
records, counts = [], []
for p in sorted((DATA / 'raw/MIXseq').glob('*.zip')):
    with zipfile.ZipFile(p) as z:
        for name in z.namelist():
            if not name.endswith('/classifications.csv') or name.startswith('__MACOSX/'): continue
            prefix = name.rsplit('/', 1)[0] + '/'
            meta = pd.read_csv(z.open(name))
            barcodes = pd.read_csv(z.open(prefix + 'barcodes.tsv'), header=None, sep='\t')[0]
            with z.open(prefix + 'matrix.mtx') as f:
                for line in f:
                    if not line.startswith(b'%'):
                        n_genes, n_cells, nnz = map(int, line.split()); break
            assert barcodes.is_unique and meta.barcode.is_unique, name
            assert n_cells == len(barcodes), name
            aligned = meta.set_index('barcode').reindex(barcodes)
            assert aligned.cell_quality.notna().all(), name
            normal = aligned[aligned.cell_quality.eq('normal') & aligned.DepMap_ID.notna()]
            c = normal.groupby('DepMap_ID').size().rename('normal_cells').reset_index()
            c['condition'] = prefix.rstrip('/'); c['archive'] = p.name
            counts.append(c)
            records.append({'archive': p.name, 'condition': prefix.rstrip('/'), 'matrix_genes': n_genes,
                'matrix_barcodes': n_cells, 'matrix_nonzero_entries': nnz, 'normal_cells_with_DepMap_ID': len(normal),
                'normal_cell_lines': normal.DepMap_ID.nunique(), 'cell_lines_ge20_normal_cells': int(c.normal_cells.ge(20).sum())})
pd.DataFrame(records).to_csv(OUT / 'MIXseq_all_archive_coverage.csv', index=False)
pd.concat(counts, ignore_index=True).to_csv(OUT / 'MIXseq_normal_cell_line_counts.csv', index=False)
result = {'downloaded_expression_zip_files': len(list((DATA / 'raw/MIXseq').glob('*.zip'))),
    'matrix_conditions': len(records), 'matrix_barcodes': sum(r['matrix_barcodes'] for r in records),
    'normal_cells_with_DepMap_ID': sum(r['normal_cells_with_DepMap_ID'] for r in records),
    'unique_normal_DepMap_IDs': pd.concat(counts).DepMap_ID.nunique(),
    'qualification': 'Dimensions and barcode metadata checked; not a full traversal or expression normalization. Author normal QC applied. Sensitivity in all_CL_features.rds derives from normalized GDSC/PRISM AUC, not a new dose-specific viability assay.'}
(OUT / 'MIXseq_actual_archive_audit.json').write_text(json.dumps(result, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))
print(json.dumps(result), flush=True)
