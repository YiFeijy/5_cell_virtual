#!/usr/bin/env python3
"""24h-specific exploratory cohort; verify original instances and L3 controls."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
DATA = ROOT / 'data'


def hash_text(x):
    return hashlib.sha256(str(x).encode()).hexdigest()


def assign(frame, key):
    counts = frame.drop_duplicates('biological_condition_group').groupby(key).size()
    loads = [0] * 5
    result = {}
    for k, size in sorted(counts.items(), key=lambda v: (-v[1], hash_text(f'20261107|{v[0]}'))):
        fold = min(range(5), key=lambda j: (loads[j], j))
        result[k] = fold
        loads[fold] += int(size)
    return frame[key].map(result)


def main():
    out = ROOT / 'data/training'
    out.mkdir(parents=True, exist_ok=True)
    v2 = DATA / 'reference'
    old = DATA / 'reference'
    c = pd.read_csv(v2 / 'hms_conditions.csv')
    counts = [{'stage': 'HMS_all_library_conditions', 'records': len(c), 'conditions': c.biological_condition_group.nunique()}]
    def stage(name):
        counts.append({'stage': name, 'records': len(c), 'conditions': c.biological_condition_group.nunique()})
    c = c[c.actual_RNA_row & c.parseable_smiles].copy()
    stage('CCLE_background_and_parseable_structure')
    c = c[c.signature_ids_24h.notna()].copy()
    stage('24h_signature_available_no_3h_requirement')
    sig = pd.read_csv(old / 'lincs_signatures.csv')
    # Select expression provenance only: no stale V1 function labels are reused.
    sig = sig[['sig_id', 'distil_id', 'distil_nsample', 'distil_cc_q75', 'tas', 'cell_id', 'pert_id', 'dose_uM', 'time_h']].copy()
    c = c.rename(columns={'signature_ids_24h': 'sig_id'})
    c = c.merge(sig, on=['sig_id', 'cell_id', 'pert_id', 'dose_uM'], validate='one_to_one')
    assert c.time_h.eq(24).all() and c.library_id.isin(['LJP005', 'LJP006']).all()
    assert c.sig_id.str.split('_').str[0].eq(c.library_id).all()
    inst = pd.read_csv(next((DATA / 'raw/LINCS_GSE70138').glob('*inst_info*')), sep='\t')
    assert inst.inst_id.is_unique
    refs = c[['sig_id', 'cell_id', 'pert_id', 'dose_uM', 'distil_id']].assign(inst_id=lambda d: d.distil_id.str.split('|')).explode('inst_id')
    refs = refs.merge(inst, on='inst_id', suffixes=('', '_instance'), validate='many_to_one')
    assert refs.cell_id.eq(refs.cell_id_instance).all() and refs.pert_id.eq(refs.pert_id_instance).all()
    assert np.isclose(refs.dose_uM, refs.pert_dose).all() and refs.pert_time.eq(24).all()
    assert refs.pert_dose_unit.str.lower().eq('um').all()
    controls = inst[inst.pert_type.eq('ctl_vehicle') & inst.pert_time.eq(24) & inst.det_plate.isin(refs.det_plate)].copy()
    assert len(controls) > 0
    z = np.load(old / 'lincs_expression.npz', allow_pickle=False)
    genes = z['gene_id'].astype(str)
    path = next((DATA / 'unpacked/LINCS_GSE70138').glob('*Level3*'))
    with h5py.File(path, 'r') as f:
        row_ids = f['0/META/COL/id'].asstr()[:]
        gene_ids = f['0/META/ROW/id'].asstr()[:]
        row_map = {s: i for i, s in enumerate(row_ids)}
        gene_map = {s: i for i, s in enumerate(gene_ids)}
        cols = np.array([gene_map[g] for g in genes])
        needed = sorted(set(refs.inst_id) | set(controls.inst_id), key=row_map.__getitem__)
        assert len(needed) == len(set(needed))
        extracted = np.empty((len(needed), len(cols)), dtype=np.float32)
        mat = f['0/DATA/0/matrix']
        for start in range(0, len(needed), 256):
            names = needed[start:start+256]
            block = mat[[row_map[s] for s in names], :]
            extracted[start:start+len(names)] = block[:, cols]
        print(json.dumps({'L3_rows_read': len(needed), 'genes': len(genes)}), flush=True)
    read_map = {s: i for i, s in enumerate(needed)}
    ctrl_info = []
    medians = {}
    for plate, group in controls.groupby('det_plate'):
        a = extracted[[read_map[s] for s in group.inst_id]]
        valid = np.isfinite(a).all(1) & (a != -666).all(1)
        good = a[valid]
        if len(good):
            medians[plate] = np.median(good, axis=0)
        ctrl_info.append({'det_plate': plate, 'listed_vehicle_instances': len(a), 'finite_non_sentinel_vehicle_instances': int(valid.sum()), 'vehicle_median_gene_sd': float(np.median(np.std(good, axis=0))) if len(good) > 1 else None})
    ctrl = pd.DataFrame(ctrl_info).set_index('det_plate')
    refs['L3_finite_non_sentinel'] = [bool(np.isfinite(extracted[read_map[s]]).all() and (extracted[read_map[s]] != -666).all()) for s in refs.inst_id]
    refs['valid_vehicle_count_same_plate'] = refs.det_plate.map(ctrl.finite_non_sentinel_vehicle_instances)
    refs['instance_technical_valid'] = refs.L3_finite_non_sentinel & refs.valid_vehicle_count_same_plate.ge(2)
    quality = []
    for sid, group in refs.groupby('sig_id'):
        good = group[group.instance_technical_valid]
        delta = np.array([extracted[read_map[r.inst_id]] - medians[r.det_plate] for r in good.itertuples()])
        correlations = []
        for i in range(len(delta)):
            for j in range(i+1, len(delta)):
                if np.std(delta[i]) > 0 and np.std(delta[j]) > 0:
                    correlations.append(float(np.corrcoef(delta[i], delta[j])[0, 1]))
        quality.append({'sig_id': sid, 'n_listed_instances': len(group), 'n_technical_valid_instances': len(good),
                        'replicate_vehicle_centered_pearson_median': float(np.median(correlations)) if correlations else np.nan,
                        'all_instances_technical_valid': bool(group.instance_technical_valid.all())})
    q = pd.DataFrame(quality)
    c = c.merge(q, on='sig_id', validate='one_to_one')
    assert np.array_equal(c.distil_nsample.to_numpy(), c.n_listed_instances.to_numpy())
    c['main_eligible'] = c.all_instances_technical_valid & c.n_technical_valid_instances.ge(2)
    c['eligibility_reason'] = np.where(c.main_eligible, 'source_and_numeric_checks_passed;at_least_two_instances_with_same_plate_controls', 'fewer_than_two_valid_instances_or_control_issue')
    c.to_csv(out / 'all_24h_candidates_with_QC.csv', index=False)
    c = c[c.main_eligible].sort_values('library_condition_id').reset_index(drop=True)
    stage('basic_technical_QC_and_at_least_two_instances')
    prior = pd.read_csv(DATA / 'reference/prior_condition_ids.csv')
    drugs = pd.read_csv(DATA / 'reference/drugs.csv')
    cells = pd.read_csv(DATA / 'reference/cells.csv')
    c = c.merge(drugs[['pert_id', 'drug_index', 'chemical_split_group', 'scaffold_split_group']], on='pert_id', validate='many_to_one')
    c = c.merge(cells[['cell_id', 'ModelID', 'cell_index']], on=['cell_id', 'ModelID'], validate='many_to_one')
    c['sample_index'] = np.arange(len(c))
    c['condition_weight'] = 1. / c.groupby('biological_condition_group').library_condition_id.transform('size')
    c['sensitivity_rep3plus'] = c.n_technical_valid_instances.ge(3)
    c['sensitivity_without_CGP60474'] = c.pert_id.ne('BRD-K79090631')
    zmap = {s: i for i, s in enumerate(z['sig_id'].astype(str))}
    response = z['expression_z'][[zmap[s] for s in c.sig_id]]
    base = pd.read_csv(v2 / 'baseline.csv', index_col=0).reindex(index=cells.ModelID, columns=genes)
    prev = np.load(DATA / 'reference/input_features.npz', allow_pickle=False)
    assert np.isfinite(response).all() and np.isfinite(base.to_numpy()).all()
    assert np.array_equal(genes, prev['gene_id'])
    np.savez_compressed(out / 'arrays.npz', expression_z24=response, baseline=base.to_numpy(dtype=np.float32),
                        drug_morgan=prev['drug_morgan'], cell_index=c.cell_index.to_numpy(), drug_index=c.drug_index.to_numpy(),
                        log10_dose_uM=np.log10(c.dose_uM.to_numpy()).astype(np.float32),
                        count72=c.relative_cell_count.to_numpy(), condition_weight=c.condition_weight.to_numpy(),
                        gene_id=genes, library_condition_id=c.library_condition_id.to_numpy(dtype=str))
    c.to_csv(out / 'samples.csv', index=False)
    refs.to_csv(out / 'instance_checks.csv', index=False)
    ctrl.reset_index().to_csv(out / 'vehicle_plate_checks.csv', index=False)
    pd.DataFrame(counts).to_csv(out / 'cohort_flow.csv', index=False)
    splits = []
    for protocol, key in [('unseen_drug', 'chemical_split_group'), ('unseen_scaffold', 'scaffold_split_group')]:
        assignment = assign(c, key)
        for fold in range(5):
            d = c[['sample_index', 'library_condition_id', 'biological_condition_group', key]].rename(columns={key: 'heldout_group'}).copy()
            d['protocol'], d['fold'] = protocol, fold
            d['role'] = np.where(assignment.eq(fold), 'test', np.where(assignment.eq((fold+1)%5), 'validation', 'train'))
            assert d.groupby('heldout_group').role.nunique().eq(1).all()
            splits.append(d)
    d = c[['sample_index', 'library_condition_id', 'biological_condition_group']].copy()
    d['heldout_group'] = c.biological_condition_group
    d['protocol'], d['fold'] = 'dose_extrapolation', 0
    d['role'] = np.where(np.isclose(c.dose_uM, 10), 'test', np.where(np.isclose(c.dose_uM, 3.33), 'validation', 'train'))
    splits.append(d)
    sp = pd.concat(splits, ignore_index=True)
    assert sp.groupby(['protocol', 'fold', 'biological_condition_group']).role.nunique().eq(1).all()
    sp.to_csv(out / 'splits.csv', index=False)
    sp.groupby(['protocol', 'fold', 'role']).agg(records=('sample_index', 'size'), conditions=('biological_condition_group', 'nunique')).reset_index().to_csv(out / 'split_counts.csv', index=False)
    sources = [v2 / 'hms_conditions.csv', old / 'lincs_signatures.csv', old / 'lincs_expression.npz', DATA / 'reference/drugs.csv', ROOT / 'docs/framework.docx']
    source_hashes = [{'path': str(p.relative_to(ROOT)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
    result = {'stage': 'exploratory_basic_technical_QC_frozen', 'expression_time_h': 24, 'function_time_h': 72,
              'records': len(c), 'distinct_conditions': c.biological_condition_group.nunique(), 'cells': c.cell_id.nunique(), 'drugs': c.pert_id.nunique(),
              'expression_shape': list(response.shape), 'new_24h_conditions_not_requiring_3h': int((~c.library_condition_id.isin(prior.library_condition_id)).sum()),
              'original_instances_checked': len(refs), 'vehicle_plates_checked': len(ctrl),
              'vehicle_instances_checked': len(controls), 'excluded_24h_candidates': int((~q.all_instances_technical_valid | q.n_technical_valid_instances.lt(2)).sum()),
              'single_instance_24h_signatures': int(q.n_listed_instances.eq(1).sum()),
              'replicate_count_distribution': c.n_technical_valid_instances.value_counts().sort_index().to_dict(),
              'sensitivity_rep3plus_records': int(c.sensitivity_rep3plus.sum()), 'sensitivity_without_CGP_records': int(c.sensitivity_without_CGP60474.sum()),
              'sum_condition_weights': float(c.condition_weight.sum()), 'no_TAS_or_response_magnitude_exclusion': True,
              'source_hashes': source_hashes,
              'qualification': 'Basic identity, numeric and replicate/control coverage eligibility; no claim of independently validated viable-cell state or complete biological quality.'}
    (out / 'audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
