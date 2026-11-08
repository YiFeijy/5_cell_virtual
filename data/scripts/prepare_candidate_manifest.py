#!/usr/bin/env python3
"""Publish availability-only candidate conditions with traceable signature IDs."""
from pathlib import Path
import json
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / 'analysis_20261106'
keys = ['cell_id', 'pert_id', 'dose_uM']
functional = pd.read_csv(OUT / 'hms_functional_conditions.csv')
signatures = pd.read_csv(OUT / 'hms_lincs_actual_signatures.csv')
baseline = pd.read_csv(OUT / 'hms_actual_baseline_coverage.csv')
structures = pd.read_csv(OUT / 'hms_phaseII_structure_coverage.csv')
manifest = functional.copy()
for time in [3, 24]:
    selected = signatures[signatures.time_h.eq(time) & signatures.finite_landmark_values]
    g = selected.groupby(keys).sig_id.agg(lambda v: '|'.join(v)).rename(f'L5_signature_ids_{time}h').reset_index()
    manifest = manifest.merge(g, on=keys, how='left', validate='one_to_one')
manifest = manifest.merge(baseline[['StrippedCellLineName','ModelID','actual_RNA_row']], left_on='cell_id', right_on='StrippedCellLineName', how='left', validate='many_to_one')
manifest = manifest.merge(structures[['pert_id','pert_iname','canonical_smiles','parseable_smiles']], on='pert_id', how='left', validate='many_to_one')
manifest['functional_time_h'] = 72
manifest['both_expression_times_available'] = manifest.L5_signature_ids_3h.notna() & manifest.L5_signature_ids_24h.notna()
manifest['inputs_and_both_times_available'] = manifest.both_expression_times_available & manifest.actual_RNA_row & manifest.parseable_smiles
manifest['pairing_scope'] = 'same_cell_drug_nominal_dose_LJP_parallel_assays; not_same_well_or_single_cell'
manifest['qualification'] = 'availability_only; expression_QC_and_structure_identity_curation_pending'
manifest.to_csv(OUT / 'HMS_3852条件_数据可用性.csv', index=False)
ready = manifest[manifest.inputs_and_both_times_available]
ready.to_csv(OUT / 'HMS_3168条件_输入与三个终点可用.csv', index=False)
r = {'available_conditions': len(ready), 'cells': ready.cell_id.nunique(), 'drugs': ready.pert_id.nunique(),
     'cell_drug_pairs': len(ready.drop_duplicates(['cell_id','pert_id'])),
     'conditions_by_cell': ready.groupby('cell_id').size().to_dict(),
     'qualification': 'Available inputs, 3h/24h expression and 72h functional endpoint; not final QC-passing training cohort.'}
(OUT / 'HMS_candidate_manifest_summary.json').write_text(json.dumps(r,indent=2))
print(json.dumps(r))
