#!/usr/bin/env python3
"""Prepare traceable HMS/LINCS relation tables, preserving endpoint identity."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[1]
OUT = DATA / 'integration_20261106'
OUT.mkdir(exist_ok=True)
OLD = DATA / 'analysis_20261106'
keys = ['cell_id','pert_id','dose_uM']
ready = pd.read_csv(OLD / 'HMS_3168条件_输入与三个终点可用.csv')
s = pd.read_csv(OLD / 'hms_lincs_actual_signatures.csv')
inst = pd.read_csv(next((DATA / 'raw/LINCS_GSE70138').glob('*inst_info*')),sep='\t')
refs = s[['sig_id','cell_id','pert_id','dose_uM','time_h','distil_id']].assign(inst_id=lambda d:d.distil_id.str.split('|')).explode('inst_id')
refs = refs.merge(inst[['inst_id','cell_id','pert_id','pert_dose','pert_dose_unit','pert_time']],on='inst_id',how='left',suffixes=('','_instance'),validate='many_to_one')
assert refs.cell_id.eq(refs.cell_id_instance).all()
assert refs.pert_id.eq(refs.pert_id_instance).all()
assert refs.pert_dose_unit.str.lower().eq('um').all()
assert np.isclose(refs.dose_uM,refs.pert_dose,rtol=1e-6,atol=1e-9).all()
assert np.isclose(refs.time_h,refs.pert_time).all()
refs.to_csv(OUT / 'HMS_LINCS_instance_condition_verification.csv',index=False)

def key_id(row):
    token=f'{row.cell_id}|{row.pert_id}|{row.dose_uM:.9g}'
    return 'HMS_' + hashlib.sha256(token.encode()).hexdigest()[:16]

ready['condition_id'] = [key_id(r) for r in ready.itertuples()]
ready['functional_condition_id'] = ready.condition_id + '_72h_cell_count'
ready['functional_readout'] = 'relative_cell_count'
ready['secondary_functional_readout'] = 'relative_net_growth'
ready = ready.rename(columns={'growth_rate':'relative_net_growth'})
ready['functional_label_origin'] = 'technical_well_mean_within_biorepeat_then_equal_biorepeat_mean'
ready['state_reference_source'] = 'CCLE_DepMap_24Q4; cross_experiment_cell_background'
ready['relation_grade'] = 'same_project_parallel_assays_candidate'
ready['split_group_cell_drug'] = ready.ModelID + '__' + ready.pert_id
ready.to_csv(OUT / 'HMS_core_condition_manifest.csv',index=False)
edges=[]
for row in ready.itertuples():
    for time in [3,24]:
        ids=getattr(row,f'L5_signature_ids_{time}h').split('|')
        edges.append({'condition_id':row.condition_id,'expression_time_h':time,
                      'signature_ids':'|'.join(ids),'n_signature_records':len(ids),
                      'functional_condition_id':row.functional_condition_id,
                      'functional_time_h':72,'functional_readout':'relative_cell_count',
                      'functional_value':row.relative_cell_count,
                      'relative_net_growth':row.relative_net_growth,
                      'relation_grade':row.relation_grade,
                      'endpoint_reference_weight_if_flattened':.5})
pd.DataFrame(edges).to_csv(OUT / 'HMS_expression_function_relations.csv',index=False)
bio = pd.read_csv(OLD / 'hms_biological_repeat_means.csv').merge(ready[keys + ['condition_id','functional_condition_id']],on=keys,how='inner',validate='many_to_one').rename(columns={'growth_rate':'relative_net_growth'})
bio['readout_origin'] = 'within_biological_repeat_technical_well_aggregate'
bio.to_csv(OUT / 'HMS_core_functional_biorepeats.csv',index=False)
result={'core_conditions':len(ready),'cells':ready.cell_id.nunique(),'drugs':ready.pert_id.nunique(),
        'expression_time_relations':len(edges),'functional_biorepeat_records':len(bio),
        'checked_instance_references':len(refs),'instance_cell_drug_dose_time_agreement':True,
        'qualification':'These are availability and provenance-check pilots; biological expression QC and chemical identity curation remain part of training-set preparation. Functional loss is intended once per condition. The .5 relation weight illustrates endpoint reuse across two times, not a biological uncertainty weight.'}
(OUT / 'HMS_relation_pilot_summary.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
