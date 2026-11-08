#!/usr/bin/env python3
"""Export a worked example from source measurements, not model predictions."""
from pathlib import Path
import json
import sys
import h5py
import numpy as np
import pandas as pd

DATA=Path(__file__).resolve().parents[1]
ROOT=DATA.parent
OUT=DATA/'integration_20261106'/'BT20_afatinib_0.04uM实例'
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'数据可行性调研_20261105/evidence/xlrd-2.0.2-py2.py3-none-any.whl'))
records=[]
for p in sorted((DATA/'raw/HMS_LINCS').glob('*.xls')):
    d=pd.read_excel(p,sheet_name='Data')
    d=d[d.CellLine.eq('BT20') & d.BRDid.eq('BRD-A58767537') & d['Small Mol Concentration'].eq(.04)].copy()
    d['source_file']=p.name
    records.append(d)
raw=pd.concat(records,ignore_index=True)
raw.to_csv(OUT/'HMS原始9个孔.csv',index=False)
bio=raw.groupby('source_file').agg(n_technical_wells=('Relative Cell Count','size'),
    relative_cell_count_mean=('Relative Cell Count','mean')).reset_index()
bio.to_csv(OUT/'HMS三批均值.csv',index=False)
sigpath=next((DATA/'raw/LINCS_GSE70138').glob('*sig_info*'))
sig=pd.read_csv(sigpath,sep='\t')
sig=sig[sig.sig_id.isin(['LJP006_BT20_3H:G24','LJP006_BT20_24H:G24'])].copy()
sig.to_csv(OUT/'LINCS两条签名说明.csv',index=False)
inst=pd.read_csv(next((DATA/'raw/LINCS_GSE70138').glob('*inst_info*')),sep='\t')
ids=[x for values in sig.distil_id for x in values.split('|')]
inst[inst.inst_id.isin(ids)].to_csv(OUT/'LINCS六条原实例说明.csv',index=False)
gi=pd.read_csv(next((DATA/'raw/LINCS_GSE70138').glob('*gene_info*')),sep='\t')
gene=gi[gi.pr_is_lm.eq(1)][['pr_gene_id','pr_gene_symbol']].copy()
gene['pr_gene_id']=gene.pr_gene_id.astype(str)
base=pd.read_csv(DATA/'analysis_20261106/hms_baseline_landmark_genes.csv',index_col=0)
gene['baseline_log2_TPM_plus_1']=[base.loc['ACH-000536',g] for g in gene.pr_gene_id]
with h5py.File(next((DATA/'unpacked/LINCS_GSE70138').glob('*Level5*.gctx'))) as f:
    gidx={g:i for i,g in enumerate(f['0/META/ROW/id'][:].astype(str))}
    sidx={s:i for i,s in enumerate(f['0/META/COL/id'][:].astype(str))}
    for time in [3,24]:
        row=f['0/DATA/0/matrix'][sidx[f'LJP006_BT20_{time}H:G24'],:]
        gene[f'expression_signature_z_{time}h']=[float(row[gidx[g]]) for g in gene.pr_gene_id]
gene.to_csv(OUT/'实际对接的978个基因.csv',index=False)
core=pd.read_csv(DATA/'integration_20261106/HMS_core_condition_manifest.csv')
one=core[core.cell_id.eq('BT20') & core.pert_id.eq('BRD-A58767537') & core.dose_uM.eq(.04)]
one.to_csv(OUT/'最终条件记录.csv',index=False)
assert len(raw)==9 and bio.n_technical_wells.eq(3).all()
assert np.isclose(bio.relative_cell_count_mean.mean(),one.relative_cell_count.iloc[0])
assert len(gene)==978 and gene.pr_gene_id.is_unique
all_conditions=pd.read_csv(DATA/'analysis_20261106/HMS_3852条件_数据可用性.csv')
cells=all_conditions.groupby('cell_id').agg(HMS_functional_conditions=('pert_id','size'),
    HMS_drugs=('pert_id','nunique'),baseline_available=('actual_RNA_row','first')).reset_index()
counts=core.groupby('cell_id').size();cells['core_complete_input_conditions']=cells.cell_id.map(counts).fillna(0).astype(int)
models=pd.read_csv(DATA/'raw/CCLE_DepMap_24Q4/Model.csv')
core_models=models[models.StrippedCellLineName.isin(cells.cell_id)]
cells=cells.merge(core_models[['StrippedCellLineName','ModelID','CellLineName','RRID']],left_on='cell_id',right_on='StrippedCellLineName',how='left',validate='one_to_one')
cells.to_csv(OUT/'HMS核心细胞系清单.csv',index=False)
cross=pd.read_csv(DATA/'integration_20261106/LINCS_CTRP_nearest_measured_candidates.csv')
cross=cross[cross.structure_relation.isin(['same_canonical_isomeric_structure','same_standardized_parent']) & cross.dose_log10_difference.le(.1)]
crosscells=cross.groupby('ModelID_candidate').agg(expression_query_candidates=('query_id','nunique'),drugs=('pert_id','nunique')).reset_index()
crosscells=crosscells.merge(models[['ModelID','CellLineName','RRID']],left_on='ModelID_candidate',right_on='ModelID',validate='one_to_one')
crosscells['identity_status']='unique_name_candidates; stable_identity_review_pending'
crosscells.to_csv(OUT/'跨库48个候选细胞系清单.csv',index=False)
summary={'cell':'BT20','ModelID':'ACH-000536','drug':'afatinib','pert_id':'BRD-A58767537','dose_uM':.04,
    'raw_HMS_well_records':len(raw),'biological_replicates':len(bio),
    'mean_relative_cell_count':float(bio.relative_cell_count_mean.mean()),
    'functional_time_h':72,'LINCS_signature_ids':sig.sig_id.tolist(),
    'gene_rows':len(gene),'scope':'same_project_parallel_condition_link; baseline from separate CCLE study; HMS and LINCS replicate sets independently preserved'}
(OUT/'实例核对.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False))
