#!/usr/bin/env python3
"""Trace a missing CCLE row and record the real completion level of integration."""
from pathlib import Path
import gzip
import hashlib
import json
import re

import h5py
import numpy as np
import pandas as pd

DATA=Path(__file__).resolve().parents[1]
OUT=DATA/'investigation_20261106'
OUT.mkdir(exist_ok=True)
BASE=DATA/'raw/CCLE_DepMap_24Q4'
release=json.loads((DATA.parent/'数据可行性调研_20261105/evidence/figshare_depmap24q4.json').read_text())
fidx={f['name']:f for f in release['files']}
checks=[]
paths=[(BASE/n,n) for n in ['Model.csv','OmicsProfiles.csv','OmicsDefaultModelProfiles.csv','OmicsExpressionProteinCodingGenesTPMLogp1.csv']]
paths += [(OUT/('DepMap24Q4_'+n),n) for n in ['README.txt','ModelCondition.csv','OmicsDefaultModelConditionProfiles.csv']]
for p,name in paths:
    h=hashlib.md5()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    expected=fidx[name].get('computed_md5') or fidx[name].get('supplied_md5')
    assert h.hexdigest()==expected,(name,h.hexdigest(),expected)
    checks.append({'file':str(p.relative_to(DATA)),'release_file':name,'bytes':p.stat().st_size,
                   'md5':h.hexdigest(),'official_md5':expected,'match':True})
(OUT/'CCLE_file_integrity_recheck.json').write_text(json.dumps(checks,indent=2))
models=pd.read_csv(BASE/'Model.csv');profiles=pd.read_csv(BASE/'OmicsProfiles.csv')
default=pd.read_csv(BASE/'OmicsDefaultModelProfiles.csv')
conditions=pd.read_csv(OUT/'DepMap24Q4_ModelCondition.csv')
exprids=pd.read_csv(BASE/'OmicsExpressionProteinCodingGenesTPMLogp1.csv',usecols=[0]).iloc[:,0]
mcf=models[models.StrippedCellLineName.eq('MCF10A')]
assert len(mcf)==1
model=mcf.iloc[0].ModelID
assert model=='ACH-001357'
mcf.to_csv(OUT/'MCF10A_DepMap_model_record.csv',index=False)
noncancer=models[models.OncotreePrimaryDisease.eq('Non-Cancerous')].copy()
noncancer['has_expression_row']=noncancer.ModelID.isin(exprids)
noncancer[['ModelID','CellLineName','has_expression_row']].to_csv(OUT/'NonCancerous_model_expression_coverage.csv',index=False)

raw=DATA/'raw/LINCS_GSE70138'
inst=pd.read_csv(next(raw.glob('*inst_info*')),sep='\t')
vehicle=inst[inst.cell_id.eq('MCF10A') & inst.pert_type.eq('ctl_vehicle') & inst.inst_id.str.startswith('LJP')].copy()
genes=pd.read_csv(next(raw.glob('*gene_info*')),sep='\t')
lm=genes[genes.pr_is_lm.eq(1)].copy();lm['pr_gene_id']=lm.pr_gene_id.astype(str)
with h5py.File(next((DATA/'unpacked/LINCS_GSE70138').glob('*Level3*.gctx'))) as f:
    sids=f['0/META/COL/id'][:].astype(str);gids=f['0/META/ROW/id'][:].astype(str)
    smap={s:i for i,s in enumerate(sids)};gmap={g:i for i,g in enumerate(gids)}
    assert vehicle.inst_id.isin(smap).all()
    vehicle['matrix_row']=vehicle.inst_id.map(smap);vehicle=vehicle.sort_values('matrix_row')
    a=f['0/DATA/0/matrix'][vehicle.matrix_row.to_numpy(),:][:,[gmap[g] for g in lm.pr_gene_id]]
    assert np.isfinite(a).all() and (a != -666).all()
np.savez_compressed(OUT/'MCF10A_LJP_vehicle_landmark_values.npz',expression=a,
                   inst_id=vehicle.inst_id.to_numpy(dtype=str),gene_id=lm.pr_gene_id.to_numpy(dtype=str))
vehicle.to_csv(OUT/'MCF10A_LJP_vehicle_instances.csv',index=False)

rnafolder=OUT/'MCF10A_RNAseq_GSE152410'
filelist=pd.read_csv(OUT/'GSE152410_filelist.txt',sep='\t')
sizes=dict(zip(filelist.Name,filelist.Size))
samplemeta=json.loads((OUT/'GSE152410_sample_metadata.json').read_text())
controls=[s for s in samplemeta if s.get('!Sample_library_strategy')==['RNA-Seq'] and any('ctrl_0' in t.lower() for t in s.get('!Sample_title',[]))]
rna_checks=[];gene_frames=[]
for sample in controls:
    sid=sample['accession'];url=sample['!Sample_supplementary_file_1'][0].replace('ftp://','https://')
    expected=int(sizes[url.rsplit('/',1)[1]])
    p=rnafolder/(sid+'_abundance.tsv.gz');part=p.with_name(p.name+'.part')
    if not p.exists() and part.exists():
        assert part.stat().st_size==expected,(sid,'partial file size incomplete')
        with gzip.open(part,'rb') as f:
            for _ in iter(lambda:f.read(1024*1024),b''):pass
        part.replace(p)
    assert p.stat().st_size==expected
    d=pd.read_csv(p,sep='\t')
    assert np.isfinite(d.tpm).all() and d.tpm.ge(0).all()
    fields=d.target_id.str.split('|',expand=True)
    d['ensembl_gene_id']=fields[1].str.replace(r'\.\d+$','',regex=True)
    d['gene_symbol']=fields[5]
    g=d.groupby(['ensembl_gene_id','gene_symbol'],dropna=False).tpm.sum().rename(sid+'_TPM').reset_index()
    gene_frames.append(g)
    rna_checks.append({'accession':sid,'title':sample['!Sample_title'][0],'source_url':url,
                       'local_file':str(p.relative_to(DATA)),'bytes':expected,'transcript_rows':len(d),
                       'gene_rows':len(g),'TPM_sum':float(d.tpm.sum()),'gzip_CRC_and_size_checked':True,
                       'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
                       'metadata_characteristics':sample['!Sample_characteristics_ch1'],
                       'inclusion_status':'independent_control_RNAseq_candidate; culture_and_ID_harmonization_pending'})
g=gene_frames[0]
for b in gene_frames[1:]:g=g.merge(b,on=['ensembl_gene_id','gene_symbol'],how='outer',validate='one_to_one')
g.to_csv(rnafolder/'MCF10A_gene_TPM_four_controls.csv',index=False)
unique=g[~g.gene_symbol.duplicated(False)]
candidate=lm[['pr_gene_id','pr_gene_symbol']].merge(unique,left_on='pr_gene_symbol',right_on='gene_symbol',how='left',validate='one_to_one')
candidate['mapping_status']=np.where(candidate.ensembl_gene_id.notna(),'exact_symbol_candidate','requires_ID_or_alias_review')
candidate.to_csv(rnafolder/'landmark_gene_symbol_mapping_candidates.csv',index=False)
(rnafolder/'file_content_checks.json').write_text(json.dumps(rna_checks,ensure_ascii=False,indent=2))

availability=pd.read_csv(DATA/'analysis_20261106/HMS_3852条件_数据可用性.csv')
sig=pd.read_csv(DATA/'analysis_20261106/hms_lincs_actual_signatures.csv')
core=pd.read_csv(DATA/'integration_20261106/HMS_core_condition_manifest.csv')
coreids=set(s for row in core.itertuples() for t in [3,24] for s in getattr(row,f'L5_signature_ids_{t}h').split('|'))
counts={'HMS_functional_conditions':len(availability),
        'conditions_both_3h_24h':int(availability.both_expression_times_available.sum()),
        'both_times_and_pinned_CCLE_baseline':int((availability.both_expression_times_available & availability.actual_RNA_row).sum()),
        'both_times_baseline_and_parseable_structure':len(core),
        'MCF10A_both_expression_times':int(availability[availability.cell_id.eq('MCF10A')].both_expression_times_available.sum()),
        'extracted_all_HMS_linked_L5_signatures':len(sig),'landmark_genes':978,
        'signatures_used_by_core_manifest':len(coreids),'core_expression_time_relations':2*len(core),
        'core_functional_biological_records':3*len(core),
        'core_final_QC_split_training':False,'cross_library_expression_matrix_extraction_complete':False}
result={'model_id':model,'RRID':mcf.iloc[0].RRID,'model_metadata_matches':len(mcf),
        'expression_row_matches':int(exprids.eq(model).sum()),
        'all_omics_profile_matches':int(profiles.ModelID.eq(model).sum()),
        'default_profile_matches':int(default.ModelID.eq(model).sum()),
        'model_condition_matches':int(conditions.ModelID.eq(model).sum()),
        'Model_csv_definition':'All models referenced by any dataset in the DepMap portal, per official 24Q4 README; not an RNA coverage table.',
        'local_download_corruption_detected':False,'wrong_primary_model_ID_detected':False,
        'RNA_release_coverage_finding':'No public RNA profile for ACH-001357 in pinned 24Q4 mapping files.',
        'underlying_unpublished_experiment_or_release_reason':'Not stated in downloaded official release README or mapping tables; cannot distinguish never generated versus other upstream choices from these files.',
        'noncancer_models_in_release':len(noncancer),'noncancer_models_with_expression':int(noncancer.has_expression_row.sum()),
        'MCF10A_LJP_vehicle_profiles':len(vehicle),'vehicle_profiles_by_time':vehicle.pert_time.value_counts().to_dict(),
        'vehicle_landmark_matrix_shape':list(a.shape),'vehicle_values_finite_and_not_missing_sentinel':True,
        'independent_RNAseq_controls_downloaded':len(rna_checks),
        'independent_gene_TPM_rows':len(g),
        'landmark_exact_symbol_candidates':int(candidate.ensembl_gene_id.notna().sum()),
        'independent_source':'GSE152410; Gross et al. 2022, doi:10.1038/s42003-022-03975-9',
        'independent_culture_note':'Different study; experimental medium omits EGF and insulin. Audit compatibility with LJP before using as background; not inserted into CCLE matrix.',
        'integration_counts_v1':counts,
        'active_core_version':'v2_库板来源核对',
        'integration_v2_summary':json.loads((DATA/'integration_20261106/v2_库板来源核对/v2实际进度.json').read_text()) if (DATA/'integration_20261106/v2_库板来源核对/v2实际进度.json').exists() else None}
(OUT/'MCF10A原因核查与对接进度.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
print(json.dumps(result,ensure_ascii=False,indent=2))
