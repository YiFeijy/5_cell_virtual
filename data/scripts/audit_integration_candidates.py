#!/usr/bin/env python3
"""Audit candidate LINCS/CTRP links using real metadata and dose observations.

Cell-name correspondences are candidates until identity curation. Drug links keep
Broad-ID evidence and independent structure comparisons. Outputs are coverage
audits, not a final QC-passing training set.
"""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

DATA = Path(__file__).resolve().parents[1]
OUT = DATA / 'integration_20261106'
OUT.mkdir(exist_ok=True)
RDLogger.DisableLog('rdApp.*')
RESULT = {'scope': 'Candidate condition links from source IDs, names, structure strings, and actual CTRP measurements; cell identity curation and expression QC pending.'}


def norm(s):
    return re.sub('[^A-Z0-9]', '', str(s).upper())


def read_lincs(source, kind):
    return pd.read_csv(next((DATA / 'raw' / source).glob('*' + kind + '*.gz')), sep='\t', low_memory=False)


def save():
    (OUT / 'candidate_audit.json').write_text(json.dumps(RESULT, ensure_ascii=False, indent=2,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))


def structure(s):
    if not isinstance(s, str) or s in ['', '-666', 'restricted']: return (None, None)
    m = Chem.MolFromSmiles(s)
    if m is None: return (None, None)
    exact = Chem.MolToSmiles(m, isomericSmiles=True)
    parent = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))
    parent = rdMolStandardize.Uncharger().uncharge(parent)
    return exact, Chem.MolToSmiles(parent, isomericSmiles=True)


def main():
    models = pd.read_csv(DATA / 'raw/CCLE_DepMap_24Q4/Model.csv')
    models['name_key'] = models.StrippedCellLineName.map(norm)
    unique_models = models[~models.name_key.duplicated(False)].copy()
    lookup = unique_models.set_index('name_key')
    ctbase = DATA / 'unpacked/CTRPv2_2015'
    cells = pd.read_csv(ctbase / 'v20.meta.per_cell_line.txt', sep='\t')
    cells['name_key'] = cells.ccl_name.map(norm)
    cells['ModelID_candidate'] = cells.name_key.map(lookup.ModelID)
    cells['mapping_evidence'] = 'unique_normalized_name_to_DepMap; RRID_vendor_curation_pending'
    cells.to_csv(OUT / 'CTRP_cell_identity_candidates.csv', index=False)
    ex = pd.read_csv(ctbase / 'v20.meta.per_experiment.txt', sep='\t')
    assert ex.groupby('experiment_id').master_ccl_id.nunique().le(1).all()
    exmap = ex[['experiment_id', 'master_ccl_id']].drop_duplicates().merge(
        cells[['master_ccl_id', 'ModelID_candidate']], on='master_ccl_id', validate='many_to_one')
    compounds = pd.read_csv(ctbase / 'v20.meta.per_compound.txt', sep='\t')
    compounds = compounds[~compounds.cpd_name.str.contains(r'\([^)]*mol/mol\)', na=False, regex=True)].copy()
    compounds['structure_exact'], compounds['structure_parent'] = zip(*compounds.cpd_smiles.map(structure))
    frames, drug_maps, cell_maps = [], [], []
    for source in ['LINCS_GSE92742', 'LINCS_GSE70138']:
        print(source, flush=True)
        c = read_lincs(source, 'cell_info')
        c['name_key'] = c.cell_id.map(norm)
        c['ModelID_candidate'] = c.name_key.map(lookup.ModelID)
        c['DepMap_RRID'] = c.name_key.map(lookup.RRID)
        c['DepMap_catalog'] = c.name_key.map(lookup.CatalogNumber)
        c['catalog_exact_agreement'] = c.provider_catalog_id.notna() & c.DepMap_catalog.notna() & c.provider_catalog_id.astype(str).eq(c.DepMap_catalog.astype(str))
        c['source'] = source
        cell_maps.append(c)
        p = read_lincs(source, 'pert_info')
        p = p[p.pert_id.isin(compounds.broad_cpd_id)].copy()
        p['lincs_structure_exact'], p['lincs_structure_parent'] = zip(*p.canonical_smiles.map(structure))
        pm = p.merge(compounds[['master_cpd_id','broad_cpd_id','cpd_name','structure_exact','structure_parent']], left_on='pert_id', right_on='broad_cpd_id', validate='one_to_one')
        pm['structure_relation'] = np.select([
            pm.lincs_structure_exact.notna() & pm.structure_exact.notna() & pm.lincs_structure_exact.eq(pm.structure_exact),
            pm.lincs_structure_parent.notna() & pm.structure_parent.notna() & pm.lincs_structure_parent.eq(pm.structure_parent),
            pm.lincs_structure_parent.isna() | pm.structure_parent.isna()],
            ['same_canonical_isomeric_structure','same_standardized_parent','structure_missing'], default='structure_conflict_review')
        pm['source'] = source
        drug_maps.append(pm)
        s = read_lincs(source, 'sig_info')
        s = s[s.pert_type.eq('trt_cp') & s.pert_id.isin(pm.pert_id)].copy()
        if 'pert_dose' in s:
            unit = s.pert_dose_unit.astype(str).str.replace('µ','u').str.replace('μ','u').str.lower()
            s['dose_uM'] = pd.to_numeric(s.pert_dose, errors='coerce') * unit.map({'um':1., 'nm':.001,'mm':1000.,'m':1e6})
            s['dose_evidence'] = 'signature_numeric_pert_dose'
        else:
            # Recover quantitative doses from every referenced instance rather than
            # treating the signature's displayed dose as a precise numeric dose.
            inst = read_lincs(source, 'inst_info')
            unit = inst.pert_dose_unit.astype(str).str.replace('µ','u').str.replace('μ','u').str.lower()
            inst['numeric_dose_uM'] = pd.to_numeric(inst.pert_dose, errors='coerce') * unit.map({'um':1.,'nm':.001,'mm':1000.,'m':1e6})
            refs = s[['sig_id','distil_id']].assign(inst_id=lambda d: d.distil_id.str.split('|')).explode('inst_id')
            refs = refs.merge(inst[['inst_id','numeric_dose_uM']], on='inst_id', how='left', validate='many_to_one')
            doses = refs.groupby('sig_id').numeric_dose_uM.agg(['first','nunique','size','count'])
            valid = doses['nunique'].eq(1) & doses['count'].eq(doses['size'])
            RESULT[source + '_dose_unresolved'] = int((~valid).sum())
            doses.loc[~valid].rename(columns={'first':'first_instance_dose_uM','nunique':'distinct_instance_doses','size':'listed_instances','count':'resolved_instances'}).to_csv(OUT / (source + '_unresolved_dose_instances.csv'))
            s['dose_uM'] = s.sig_id.map(doses.loc[valid,'first'])
            s['dose_evidence'] = 'all_distil_instances_same_numeric_dose'
        s['time_h'] = pd.to_numeric(s.pert_itime.str.extract(r'^([\d.]+)')[0], errors='coerce')
        s = s[s.dose_uM.gt(0) & s.time_h.isin([3.,6.,24.])].copy()
        s = s.merge(c[['cell_id','ModelID_candidate','catalog_exact_agreement']], on='cell_id', how='left', validate='many_to_one')
        s = s.merge(pm[['pert_id','master_cpd_id','structure_relation']],on='pert_id',validate='many_to_one')
        s['source'] = source
        frames.append(s)
    pd.concat(cell_maps).to_csv(OUT / 'LINCS_cell_identity_candidates.csv', index=False)
    drugs = pd.concat(drug_maps, ignore_index=True)
    drugs.to_csv(OUT / 'LINCS_CTRP_drug_identity_candidates.csv', index=False)
    RESULT['drug_correspondences_by_release'] = drugs.groupby(['source','structure_relation']).size().rename('rows').reset_index().to_dict('records')
    queries = pd.concat(frames, ignore_index=True).dropna(subset=['ModelID_candidate'])
    queries = queries[queries.ModelID_candidate.isin(exmap.ModelID_candidate.dropna())]
    qkeys = ['source','ModelID_candidate','pert_id','master_cpd_id','dose_uM','time_h','structure_relation']
    q = queries.groupby(qkeys, dropna=False).agg(sig_ids=('sig_id', lambda x:'|'.join(x)), n_signature_records=('sig_id','size')).reset_index()
    q['query_id'] = ['Q' + str(i) for i in range(len(q))]
    q.to_csv(OUT / 'LINCS_expression_condition_candidates.csv', index=False)
    RESULT['unique_expression_queries_before_actual_functional_coverage'] = len(q)
    print('Reading actual CTRP dose rows', flush=True)
    expids = set(exmap.loc[exmap.ModelID_candidate.isin(q.ModelID_candidate),'experiment_id'])
    compoundids = set(q.master_cpd_id)
    parts = []
    for d in pd.read_csv(ctbase / 'v20.data.per_cpd_post_qc.txt',sep='\t',chunksize=500000,
                         usecols=['experiment_id','master_cpd_id','cpd_conc_umol','cpd_avg_pv']):
        d = d[d.experiment_id.isin(expids) & d.master_cpd_id.isin(compoundids)]
        if len(d): parts.append(d)
    f = pd.concat(parts,ignore_index=True).merge(exmap,on='experiment_id',validate='many_to_one')
    f = f[np.isfinite(f.cpd_avg_pv) & f.cpd_conc_umol.gt(0)].copy()
    assert not f.duplicated(['experiment_id','master_cpd_id','cpd_conc_umol']).any()
    grouped_q = dict(tuple(q.groupby(['ModelID_candidate','master_cpd_id'])))
    records=[]
    for (model,compound,experiment), curve in f.groupby(['ModelID_candidate','master_cpd_id','experiment_id']):
        sub=grouped_q.get((model,compound))
        if sub is None: continue
        curve=curve.sort_values('cpd_conc_umol')
        ds=curve.cpd_conc_umol.to_numpy();vs=curve.cpd_avg_pv.to_numpy(); logs=np.log10(ds)
        for row in sub.itertuples():
            j=int(np.abs(logs-np.log10(row.dose_uM)).argmin())
            records.append({'query_id':row.query_id,'experiment_id':experiment,
                'CTRP_curve_id':f'CTRP2015:{experiment}:{compound}',
                'functional_dose_uM':ds[j], 'measured_relative_viability':vs[j],
                'dose_log10_difference':abs(float(logs[j]-np.log10(row.dose_uM))),
                'numeric_dose_match':bool(np.isclose(ds[j],row.dose_uM,rtol=1e-6,atol=1e-9)),
                'within_measured_dose_range':bool(ds[0]<=row.dose_uM<=ds[-1] and len(ds)>=2),
                'curve_min_dose_uM':ds[0],'curve_max_dose_uM':ds[-1],
                'functional_time_h':72,'functional_time_evidence':'CTRP2015_protocol; no per-row duration field',
                'label_origin':'replicate_aggregate','cell_identity_status':'candidate_name_mapping'})
    links=pd.DataFrame(records).merge(q,on='query_id',validate='many_to_one')
    links.to_csv(OUT / 'LINCS_CTRP_nearest_measured_candidates.csv',index=False)
    counts=[]
    for identity,label in [(links,'Broad_ID_candidates'),(links[links.structure_relation.isin(['same_canonical_isomeric_structure','same_standardized_parent'])],'Broad_ID_and_structure_supported')]:
        for scope,s in [('all_releases',identity), *list(identity.groupby('source'))]:
            for mode,mask in [('numeric_dose_match',s.numeric_dose_match),('log10_delta_le_0.05',s.dose_log10_difference.le(.05)),('log10_delta_le_0.1',s.dose_log10_difference.le(.1)),('log10_delta_le_0.2',s.dose_log10_difference.le(.2)),('within_measured_range',s.within_measured_dose_range)]:
                a=s[mask];counts.append({'identity_scope':label,'release_scope':scope,'dose_rule':mode,
                    'relation_rows':len(a),'unique_expression_queries':a.query_id.nunique(),
                    'cell_models':a.ModelID_candidate.nunique(),'drugs':a.pert_id.nunique(),
                    'functional_curves':a.CTRP_curve_id.nunique(),
                    'cell_drug_pairs':len(a.drop_duplicates(['ModelID_candidate','pert_id']))})
    pd.DataFrame(counts).to_csv(OUT / 'LINCS_CTRP_coverage_by_matching_rule.csv',index=False)
    RESULT['coverage'] = counts
    RESULT['qualification'] = 'All counts are condition-link candidates, not independent sample counts or final training eligibility. Same Broad ID plus compatible structure still requires cell-identity review. Matching dose across studies preserves separate assay and time; no expression-derived or interpolated viability labels generated.'
    import rdkit
    RESULT['rdkit_version'] = rdkit.__version__
    RESULT['structure_standardization'] = 'Canonical isomeric SMILES; comparison additionally uses Cleanup, FragmentParent and Uncharger; no stereo removal or tautomer merging.'
    save()
    print(json.dumps(RESULT,ensure_ascii=False),flush=True)


if __name__ == '__main__': main()
