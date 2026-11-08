#!/usr/bin/env python3
"""V2 candidate preparation: use the actual HMS library key before aggregation.

Keeps V1 as a broader condition-candidate audit. No training or biological-QC
eligibility is implied by these preparation tables.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

DATA=Path(__file__).resolve().parents[1]
ROOT=DATA.parent
OUT=DATA/'integration_20261106/v2_库板来源核对'
OUT.mkdir(exist_ok=True)
OLD=DATA/'analysis_20261106'
sys.path.insert(0,str(ROOT/'数据可行性调研_20261105/evidence/xlrd-2.0.2-py2.py3-none-any.whl'))
frames=[]
for p in sorted((DATA/'raw/HMS_LINCS').glob('*.xls')):
    d=pd.read_excel(p,sheet_name='Data')
    d['source_excel_row']=np.arange(len(d))+2
    d['source_file']=p.name
    frames.append(d)
w=pd.concat(frames,ignore_index=True).rename(columns={'CellLine':'cell_id','BRDid':'pert_id',
    'Small Mol Concentration':'dose_uM','Relative Cell Count':'relative_cell_count',
    'Relative Growth Rate After Treatment':'relative_net_growth'})
w['library_id']=w['LJP Library Plate'].astype(str).str.extract(r'(\d+)')[0].astype(int).map(lambda x:f'LJP{x:03d}')
w['dose_uM']=pd.to_numeric(w.dose_uM)
assert np.isfinite(w[['relative_cell_count','relative_net_growth','dose_uM']].to_numpy()).all()
assert w['Small Mol Conc Unit'].str.lower().eq('um').all() and w.Timepoint.eq(72).all()
keys=['library_id','cell_id','pert_id','dose_uM']
w['well_reference']=w.source_file+':'+w.source_excel_row.astype(str)
bio=w.groupby(keys+['source_file']).agg(relative_cell_count=('relative_cell_count','mean'),
    relative_net_growth=('relative_net_growth','mean'),n_technical_rows=('well_reference','size'),
    original_well_references=('well_reference',lambda x:'|'.join(x))).reset_index()
fun=bio.groupby(keys).agg(relative_cell_count=('relative_cell_count','mean'),
    relative_net_growth=('relative_net_growth','mean'),n_bio=('source_file','nunique'),
    relative_cell_count_sd=('relative_cell_count','std')).reset_index()
assert fun.n_bio.eq(3).all()

def key_id(r,library=True):
    k=f'{r.cell_id}|{r.pert_id}|{r.dose_uM:.9g}'
    if library:k=r.library_id+'|'+k
    return ('HMSv2_LIB_' if library else 'HMSv2_COND_')+hashlib.sha256(k.encode()).hexdigest()[:16]

fun['library_condition_id']=[key_id(r) for r in fun.itertuples()]
fun['biological_condition_group']=[key_id(r,False) for r in fun.itertuples()]
fun['functional_time_h']=72
fun['readout']='relative_cell_count'
fun['label_origin']='technical_mean_within_library_biorepeat_then_equal_biorepeat_mean'
sig=pd.read_csv(OLD/'hms_lincs_actual_signatures.csv')
sig=sig.drop(columns=['relative_cell_count','growth_rate','n_bio','relative_cell_count_sd'],errors='ignore')
sig['library_id']=sig.sig_id.str.split('_').str[0]
allowed=set(fun.library_id)
extra=sig[~sig.library_id.isin(allowed)].copy()
extra['reason_for_separate_pool']='library_not_in_HMS_source_plates; retained_as_broader_LJP_candidates'
extra.to_csv(OUT/'其他LJP库板签名_另行保留.csv',index=False)
strict=sig[sig.library_id.isin(allowed)].copy()
matching=strict.merge(fun[keys+['library_condition_id','biological_condition_group']],on=keys,how='left',validate='many_to_one')
matching[matching.library_condition_id.isna()].to_csv(OUT/'同库板条件仍未对应的签名.csv',index=False)
matching=matching.dropna(subset=['library_condition_id'])
assert matching.finite_landmark_values.all()
matched_times=matching.groupby('library_condition_id').time_h.agg(set)
fun['both_expression_times']=fun.library_condition_id.map(matched_times).apply(lambda x:isinstance(x,set) and {3.,24.}.issubset(x))
base=pd.read_csv(OLD/'hms_actual_baseline_coverage.csv').rename(columns={'StrippedCellLineName':'cell_id'})
struct=pd.read_csv(OLD/'hms_phaseII_structure_coverage.csv')
fun=fun.merge(base[['cell_id','ModelID','actual_RNA_row']],on='cell_id',how='left',validate='many_to_one')
fun=fun.merge(struct[['pert_id','pert_iname','canonical_smiles','parseable_smiles']],on='pert_id',how='left',validate='many_to_one')
fun['input_and_observation_available']=fun.both_expression_times & fun.actual_RNA_row & fun.parseable_smiles
fun['QC_stage']='source_library_and_numeric_condition_checked; biological_expression_QC_pending'
for t in [3,24]:
    ids=matching[matching.time_h.eq(t)].groupby('library_condition_id').sig_id.agg(lambda x:'|'.join(x))
    fun[f'signature_ids_{t}h']=fun.library_condition_id.map(ids)
core=fun[fun.input_and_observation_available].copy()
core['split_group_cell_drug']=core.ModelID+'__'+core.pert_id
fun.to_csv(OUT/'全部HMS库板条件及可用性.csv',index=False)
core.to_csv(OUT/'核心库板条件表.csv',index=False)
bio=bio.merge(fun[keys+['library_condition_id','biological_condition_group']],on=keys,validate='many_to_one')
bio.to_csv(OUT/'库板内生物重复与技术孔引用.csv',index=False)
core_sig=matching[matching.library_condition_id.isin(core.library_condition_id)].copy()
z=np.load(OLD/'hms_lincs_landmark_expression.npz')
oldidx={s:i for i,s in enumerate(z['sig_id'])}
assert core_sig.sig_id.is_unique and core_sig.sig_id.isin(oldidx).all()
core_sig=core_sig.sort_values(['library_condition_id','time_h','sig_id']).reset_index(drop=True)
core_sig['export_matrix_row']=np.arange(len(core_sig))
core_sig.to_csv(OUT/'核心表达签名与矩阵行号.csv',index=False)
arr=z['expression_z'][[oldidx[s] for s in core_sig.sig_id],:]
assert np.isfinite(arr).all() and arr.shape[1]==978
np.savez_compressed(OUT/'核心表达签名矩阵.npz',expression_z=arr,
    sig_id=core_sig.sig_id.to_numpy(dtype=str),gene_id=z['gene_id'])
baseline=pd.read_csv(OLD/'hms_baseline_landmark_genes.csv',index_col=0)
baseline=baseline.loc[baseline.index.isin(core.ModelID)].reindex(columns=z['gene_id'])
assert np.isfinite(baseline.to_numpy()).all()
baseline.to_csv(OUT/'CCLE背景矩阵_与表达基因同顺序.csv')
edge=[]
for r in core.itertuples():
    for t in [3,24]:
        edge.append({'library_condition_id':r.library_condition_id,
            'biological_condition_group':r.biological_condition_group,'library_id':r.library_id,
            'cell_id':r.cell_id,'pert_id':r.pert_id,'dose_uM':r.dose_uM,
            'expression_time_h':t,'signature_ids':getattr(r,f'signature_ids_{t}h'),
            'functional_time_h':72,'relative_cell_count':r.relative_cell_count,
            'relative_net_growth':r.relative_net_growth,'endpoint_weight_if_flattened':.5})
edges=pd.DataFrame(edge)
assert edges.groupby('library_condition_id').endpoint_weight_if_flattened.sum().eq(1).all()
edges.to_csv(OUT/'核心表达功能关系表.csv',index=False)
core.drop_duplicates('biological_condition_group').groupby('cell_id').size().rename('unique_cell_drug_dose_conditions').to_csv(OUT/'核心细胞条件数.csv')
summary={'V1_signature_pool':len(sig),'HMS_source_libraries':sorted(allowed),
    'other_library_signatures_kept_separately':len(extra),
    'other_libraries':extra.library_id.value_counts().to_dict(),
    'same_library_signature_pool':len(strict),'same_library_and_condition_matched_signatures':len(matching),
    'HMS_raw_technical_rows':len(w),'HMS_library_conditions':len(fun),
    'HMS_distinct_cell_drug_dose_conditions':fun.biological_condition_group.nunique(),
    'all_library_conditions_with_both_times':int(fun.both_expression_times.sum()),
    'core_library_condition_records':len(core),'core_distinct_cell_drug_dose_conditions':core.biological_condition_group.nunique(),
    'core_cell_models':core.cell_id.nunique(),'core_drugs':core.pert_id.nunique(),
    'core_expression_signatures':len(core_sig),'core_expression_matrix_shape':list(arr.shape),
    'core_expression_function_relation_rows':len(edges),'baseline_matrix_shape':list(baseline.shape),
    'core_functional_biorepeat_records':int(bio.library_condition_id.isin(core.library_condition_id).sum()),
    'technical_rows_per_library_biorepeat_distribution':bio.n_technical_rows.value_counts().to_dict(),
    'final_training_set_released':False,'models_trained':0,
    'interpretation':'Library records preserve repeat/design structure. Distinct biological condition counts group cell, drug and dose across libraries. Matrix extraction and source checks completed; biological-QC eligibility and split freezing remain.'}
(OUT/'v2实际进度.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2))
