#!/usr/bin/env python3
"""Inspect real downloaded measurements; descriptive audit, not model training."""
from pathlib import Path
from datetime import datetime
import collections
import json
import re
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data'
OUT = DATA / 'analysis_20261106'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / '数据可行性调研_20261105/evidence/xlrd-2.0.2-py2.py3-none-any.whl'))
RESULT = {'created_at': datetime.now().astimezone().isoformat(),
          'scope': 'Actual matrix and table inspection. Descriptive correlations, no trained model; no claim of same-well pairing.'}


def save():
    (OUT / 'audit_results.json').write_text(json.dumps(RESULT, ensure_ascii=False, indent=2,
        default=lambda v: v.item() if isinstance(v, np.generic) else str(v)))


def corr(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    return {'n': int(ok.sum()), 'spearman': float(spearmanr(a[ok], b[ok]).statistic) if ok.sum() > 2 else None}


def distribution(v):
    a = np.asarray(v, dtype=float)
    a = a[np.isfinite(a)]
    return dict(zip(['min', 'q01', 'q25', 'median', 'q75', 'q99', 'max'],
                    np.quantile(a, [0, .01, .25, .5, .75, .99, 1]).tolist())) if len(a) else {}


def read_meta(stem):
    return pd.read_csv(next((DATA / 'raw/LINCS_GSE70138').glob('*' + stem + '*.gz')), sep='\t', low_memory=False)


def functional():
    frames = []
    for p in sorted((DATA / 'raw/HMS_LINCS').glob('*.xls')):
        d = pd.read_excel(p, sheet_name='Data')
        d['bio_file'] = p.name
        frames.append(d)
    d = pd.concat(frames, ignore_index=True).rename(columns={
        'CellLine': 'cell_id', 'BRDid': 'pert_id', 'Small Mol Concentration': 'dose_uM',
        'Relative Cell Count': 'relative_cell_count',
        'Relative Growth Rate After Treatment': 'growth_rate'})
    d['dose_uM'] = pd.to_numeric(d.dose_uM).round(6)
    keys = ['cell_id', 'pert_id', 'dose_uM']
    numeric = ['relative_cell_count', 'growth_rate']
    bio = d.groupby(keys + ['bio_file'], dropna=False)[numeric].mean().reset_index()
    bio.to_csv(OUT / 'hms_biological_repeat_means.csv', index=False)
    cond = bio.groupby(keys).agg(relative_cell_count=('relative_cell_count', 'mean'),
        growth_rate=('growth_rate', 'mean'), n_bio=('bio_file', 'nunique'),
        relative_cell_count_sd=('relative_cell_count', 'std')).reset_index()
    cond.to_csv(OUT / 'hms_functional_conditions.csv', index=False)
    wide = bio.pivot(index=keys, columns='bio_file', values='relative_cell_count')
    rep_corr = []
    for i, a in enumerate(wide.columns):
        for b in wide.columns[i + 1:]:
            rep_corr.append({'a': a, 'b': b, **corr(wide[a], wide[b])})
    # Confirm that stored relative count agrees with the accompanying measured counts.
    recompute = pd.to_numeric(d['Total Cell Count After Treatment']) / pd.to_numeric(d['Total Control Cell Count'])
    RESULT['HMS'] = {'well_rows': len(d), 'conditions': len(cond), 'cells': cond.cell_id.nunique(),
        'drugs': cond.pert_id.nunique(), 'doses_uM': sorted(cond.dose_uM.unique()),
        'endpoint_h': sorted(d.Timepoint.unique()), 'bio_repeat_count_distribution': cond.n_bio.value_counts().to_dict(),
        'relative_cell_count': distribution(cond.relative_cell_count), 'growth_rate': distribution(cond.growth_rate),
        'conditions_relative_count_below_half': int(cond.relative_cell_count.lt(.5).sum()),
        'conditions_negative_growth_rate': int(cond.growth_rate.lt(0).sum()),
        'repeat_correlations': rep_corr,
        'stored_ratio_max_absolute_difference': float(np.nanmax(np.abs(recompute - d.relative_cell_count))),
        'qualification': 'Relative cell count versus control and growth rate are distinct endpoints; neither is a direct fraction of dead cells. Average technical wells within each biological file, then weight biological files equally.'}
    return cond, keys


def baseline_and_structures(cond):
    gi = read_meta('gene_info')
    landmark_ids = set(gi.loc[gi.pr_is_lm.eq(1), 'pr_gene_id'].astype(str))
    models = pd.read_csv(DATA / 'raw/CCLE_DepMap_24Q4/Model.csv')
    models = models[models.StrippedCellLineName.isin(cond.cell_id.unique())]
    exp_path = DATA / 'raw/CCLE_DepMap_24Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv'
    header = pd.read_csv(exp_path, nrows=0).columns
    gene_columns = {c: re.search(r'\((\d+)\)$', c).group(1) for c in header[1:] if re.search(r'\((\d+)\)$', c)}
    use = [header[0]] + [c for c, gid in gene_columns.items() if gid in landmark_ids]
    exp = pd.read_csv(exp_path, usecols=use, index_col=0)
    mapped = models[['ModelID', 'StrippedCellLineName']].copy()
    mapped['actual_RNA_row'] = mapped.ModelID.isin(exp.index)
    mapped.to_csv(OUT / 'hms_actual_baseline_coverage.csv', index=False)
    sub = exp.loc[exp.index.isin(mapped.ModelID)].copy()
    sub.columns = [gene_columns[c] for c in sub.columns]
    sub.to_csv(OUT / 'hms_baseline_landmark_genes.csv')
    pert = read_meta('pert_info')
    compounds = cond[['pert_id']].drop_duplicates().merge(pert, on='pert_id', how='left', validate='one_to_one')
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.error')
    compounds['parseable_smiles'] = compounds.canonical_smiles.apply(
        lambda s: isinstance(s, str) and s not in ['', '-666'] and Chem.MolFromSmiles(s) is not None)
    compounds.to_csv(OUT / 'hms_phaseII_structure_coverage.csv', index=False)
    RESULT['inputs'] = {'baseline_actual_rows': len(exp), 'baseline_gene_columns': len(header) - 1,
        'lincs_landmark_genes': len(landmark_ids), 'baseline_landmark_overlap': len(use) - 1,
        'missing_landmark_gene_ids': sorted(landmark_ids - set(gene_columns.values())),
        'hms_cells_with_actual_baseline': mapped[mapped.actual_RNA_row].StrippedCellLineName.tolist(),
        'hms_cells_without_actual_baseline': mapped[~mapped.actual_RNA_row].StrippedCellLineName.tolist(),
        'finite_baseline_values': bool(np.isfinite(sub.to_numpy()).all()),
        'hms_drugs_parseable_phaseII_smiles': int(compounds.parseable_smiles.sum()),
        'hms_drugs_missing_or_invalid_phaseII_smiles': compounds[~compounds.parseable_smiles].pert_id.tolist()}
    return landmark_ids, set(mapped.loc[mapped.actual_RNA_row, 'StrippedCellLineName']), set(compounds.loc[compounds.parseable_smiles, 'pert_id'])


def expression(cond, keys, landmark_ids, baseline_cells, structure_drugs):
    sig = read_meta('sig_info')
    metrics = read_meta('sig_metrics')
    metric_cols = ['tas', 'distil_cc_q75', 'distil_nsample', 'distil_ss']
    metrics[metric_cols] = metrics[metric_cols].replace(-666, np.nan)
    sig = sig[sig.sig_id.str.startswith('LJP') & sig.pert_type.eq('trt_cp')].copy()
    sig['dose_uM'] = pd.to_numeric(sig.pert_idose.str.extract(r'^([\d.eE+-]+)')[0], errors='coerce').round(6)
    units = sig.pert_idose.str.extract(r'\s+(\S+)')[0].str.replace('µ', 'u').str.replace('μ', 'u')
    assert units.dropna().isin(['um', 'uM']).all(), units.unique()
    sig['time_h'] = pd.to_numeric(sig.pert_itime.str.extract(r'^([\d.]+)')[0])
    matched = sig.merge(cond, on=keys, how='inner', validate='many_to_one').merge(metrics, on='sig_id', how='left', suffixes=('', '_metric'), validate='one_to_one')
    l5 = next((DATA / 'unpacked/LINCS_GSE70138').glob('*Level5*.gctx'))
    with h5py.File(l5) as f:
        gids = f['0/META/ROW/id'][:].astype(str)
        ids = f['0/META/COL/id'][:].astype(str)
        rowmap = {v: i for i, v in enumerate(ids)}
        matched['matrix_row'] = matched.sig_id.map(rowmap)
        missing_ids = int(matched.matrix_row.isna().sum())
        matched = matched.dropna(subset=['matrix_row']).sort_values('matrix_row').reset_index(drop=True)
        idx = np.flatnonzero(np.isin(gids, list(landmark_ids)))
        arr = np.empty((len(matched), len(idx)), dtype=np.float32)
        m = f['0/DATA/0/matrix']
        for start in range(0, len(matched), 256):
            rows = matched.matrix_row.iloc[start:start + 256].astype(int).to_numpy()
            arr[start:start + len(rows)] = m[rows, :][:, idx]
    finite = np.isfinite(arr).all(axis=1) & (arr != -666).all(axis=1)
    matched['finite_landmark_values'] = finite
    matched['landmark_rms_z'] = np.sqrt(np.mean(arr ** 2, axis=1))
    matched['landmark_abs_z_ge2_count'] = (np.abs(arr) >= 2).sum(axis=1)
    matched.to_csv(OUT / 'hms_lincs_actual_signatures.csv', index=False)
    np.savez_compressed(OUT / 'hms_lincs_landmark_expression.npz', expression_z=arr,
                        sig_id=matched.sig_id.to_numpy(dtype=str), gene_id=gids[idx])
    measures = ['landmark_rms_z', 'landmark_abs_z_ge2_count', 'tas', 'distil_cc_q75']
    obs = matched[finite].groupby(keys + ['time_h'])[measures].mean().reset_index().merge(cond, on=keys, validate='many_to_one')
    obs.to_csv(OUT / 'hms_expression_functional_conditions.csv', index=False)
    times = obs.groupby(keys).time_h.agg(set)
    actual_both = times[times.apply(lambda x: {3., 24.}.issubset(x))].index
    ready = obs[obs.cell_id.isin(baseline_cells) & obs.pert_id.isin(structure_drugs)]
    ready_times = ready.groupby(keys).time_h.agg(set)
    records = []
    for time, sub in obs.groupby('time_h'):
        for endpoint in ['relative_cell_count', 'growth_rate']:
            records.append({'time_h': time, 'group': 'all', 'endpoint': endpoint,
                            **corr(sub.landmark_rms_z, sub[endpoint])})
            for cell, c in sub.groupby('cell_id'):
                records.append({'time_h': time, 'group': cell, 'endpoint': endpoint,
                                **corr(c.landmark_rms_z, c[endpoint])})
    pd.DataFrame(records).to_csv(OUT / 'response_function_correlation.csv', index=False)
    wide = obs.pivot(index=keys, columns='time_h', values='landmark_rms_z').dropna()
    common = wide.merge(cond.set_index(keys)[['relative_cell_count']], left_index=True, right_index=True)
    same_cohort = {'same_condition_both_times_n': len(common),
        '3h_spearman': corr(common[3.], common.relative_cell_count)['spearman'],
        '24h_spearman': corr(common[24.], common.relative_cell_count)['spearman']}
    (OUT / 'response_function_same_cohort.json').write_text(json.dumps(same_cohort, indent=2))
    by_dose = [{'time_h': t, 'dose_uM': dose, **corr(s.landmark_rms_z, s.relative_cell_count)}
               for (t, dose), s in obs.groupby(['time_h', 'dose_uM'])]
    pd.DataFrame(by_dose).to_csv(OUT / 'response_function_by_dose.csv', index=False)
    # Every distil_id is checked against the downloaded L3 matrix and its own vehicle plate.
    inst = read_meta('inst_info').set_index('inst_id', drop=False)
    l3 = next((DATA / 'unpacked/LINCS_GSE70138').glob('*Level3*.gctx'))
    with h5py.File(l3) as f:
        l3ids = set(f['0/META/COL/id'][:].astype(str))
    vehicles = inst[inst.pert_type.eq('ctl_vehicle')]
    vehicle_plates = set(vehicles.det_plate)
    detail = []
    for row in matched.itertuples():
        srcids = str(row.distil_id).split('|')
        available = [s for s in srcids if s in l3ids]
        with_control = [s for s in available if s in inst.index and inst.at[s, 'det_plate'] in vehicle_plates]
        detail.append({'sig_id': row.sig_id, 'listed_instances': len(srcids),
                       'instances_in_L3': len(available), 'instances_with_same_plate_vehicle_metadata': len(with_control)})
    detail = pd.DataFrame(detail)
    detail.to_csv(OUT / 'hms_L3_instance_control_coverage.csv', index=False)
    RESULT['joint_actual_matrices'] = {'signature_records': len(matched), 'signature_ids_missing_from_L5': missing_ids,
        'landmark_columns_read': arr.shape[1], 'signatures_with_all_landmark_values_finite': int(finite.sum()),
        'unique_conditions_by_time': obs.groupby('time_h').size().to_dict(),
        'conditions_with_both_3h_and_24h': len(actual_both),
        'both_times_with_baseline_and_parseable_structure': int(ready_times.apply(lambda x: {3., 24.}.issubset(x)).sum()),
        'signature_TAS': distribution(matched.tas), 'signature_distil_cc_q75': distribution(matched.distil_cc_q75),
        'signature_missing_distil_cc_q75': int(matched.distil_cc_q75.isna().sum()),
        'signature_missing_TAS': int(matched.tas.isna().sum()),
        'signature_TAS_below_0_2': int(matched.tas.lt(.2).sum()),
        'signature_distil_cc_q75_below_0_2': int(matched.distil_cc_q75.lt(.2).sum()),
        'L3_total_listed_replicate_references': int(detail.listed_instances.sum()),
        'L3_total_available_replicate_references': int(detail.instances_in_L3.sum()),
        'L3_signature_all_references_available': int(detail.instances_in_L3.eq(detail.listed_instances).sum()),
        'L3_signature_all_references_have_same_plate_vehicle_metadata': int(detail.instances_with_same_plate_vehicle_metadata.eq(detail.listed_instances).sum()),
        'correlations': records,
        'qualification': 'Condition correspondence in coordinated project; independent time points and wells. Matrix presence and finite values do not prove biological quality. TAS/CC thresholds are descriptive counts, not exclusion rules. L3 control presence verified in metadata, not yet control quality or normalization. RMS z is expression response magnitude, not viability.'}
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for ax, time in zip(axes, [3., 24.]):
            q = obs[obs.time_h.eq(time)]
            for cell, c in q.groupby('cell_id'):
                ax.scatter(c.landmark_rms_z, c.relative_cell_count, s=4, alpha=.35, label=cell)
            ax.set_title(f'{int(time)} h expression / 72 h relative count')
            ax.set_xlabel('RMS z-score across 978 measured genes')
            ax.set_ylim(-.05, min(2.5, obs.relative_cell_count.max() + .05))
        axes[0].set_ylabel('72 h relative cell count (treated / control)')
        axes[1].legend(markerscale=3, fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / 'expression_vs_relative_cell_count.png', dpi=180)
        plt.close(fig)
    except ImportError:
        pass


def sciplex():
    result = {}
    feature_audit = {}
    for p in (DATA / 'raw/sciPlex').glob('*.h5ad'):
        with h5py.File(p) as f:
            def column(k):
                o = f['obs'][k]
                if isinstance(o, h5py.Dataset):
                    a = o[:]
                    return a.astype(str) if a.dtype.kind in 'SO' else a
                categories = o['categories'][:]
                codes = o['codes'][:]
                if categories.dtype.kind in 'SO': categories = categories.astype(str)
                out = np.full(len(codes), None, dtype=object)
                ok = codes >= 0
                out[ok] = categories[codes[ok]]
                return out
            keep = [k for k in ['cell_line', 'perturbation', 'dose_value', 'dose_unit', 'time', 'replicate', 'ncounts', 'nperts'] if k in f['obs']]
            obs = pd.DataFrame({k: column(k) for k in keep})
            keys = [k for k in keep if k not in ['ncounts', 'nperts']]
            design = obs.groupby(keys, dropna=False).size().rename('raw_cell_count').reset_index()
            design.to_csv(OUT / (p.stem + '_condition_raw_cell_counts.csv'), index=False)
            result[p.name] = {'X_shape': list(f['X'].attrs.get('shape', f['X'].shape if isinstance(f['X'], h5py.Dataset) else [])),
                'obs_rows': len(obs), 'cell_lines': obs.cell_line.value_counts(dropna=False).to_dict(),
                'perturbation_names': obs.perturbation.nunique(), 'conditions_with_replicate': len(design),
                'condition_cell_count': distribution(design.raw_cell_count),
                'raw_conditions_fewer_than_20_cells': int(design.raw_cell_count.lt(20).sum()),
                'missing_cell_line': int(obs.cell_line.isna().sum()), 'missing_perturbation': int(obs.perturbation.isna().sum()),
                'time_values': obs.time.value_counts(dropna=False).to_dict() if 'time' in obs else None,
                'nperts_values': obs.nperts.value_counts(dropna=False).to_dict() if 'nperts' in obs else None,
                'qualification': 'Raw file rows and annotations; cell QC and doublet assessment not performed. Captured cell numbers are not an independent viability assay.'}
            gene = f['var/ensembl_id']
            categories = gene['categories'][:].astype(str)
            codes = gene['codes'][:]
            values = categories[codes[codes >= 0]]
            feature_audit[p.name] = {'var_rows': len(codes), 'unique_ensembl_ids': len(set(values)),
                'human_prefix_ENSG_rows': int(np.char.startswith(values, 'ENSG').sum()),
                'mouse_prefix_ENSMUSG_rows': int(np.char.startswith(values, 'ENSMUSG').sum()),
                'missing_ensembl_rows': int((codes < 0).sum())}
    (OUT / 'sciPlex_feature_identity_audit.json').write_text(json.dumps(feature_audit, indent=2))
    RESULT['sciPlex'] = result


def ctrp():
    base = DATA / 'unpacked/CTRPv2_2015'
    ex = pd.read_csv(base / 'v20.meta.per_experiment.txt', sep='\t')
    compounds = pd.read_csv(base / 'v20.meta.per_compound.txt', sep='\t')
    cells = pd.read_csv(base / 'v20.meta.per_cell_line.txt', sep='\t')
    combo = compounds.cpd_name.str.contains(r'\([^)]*mol/mol\)', regex=True, na=False)
    combination_ids = set(compounds.loc[combo, 'master_cpd_id'])
    compounds.assign(explicit_fixed_ratio_combination=combo).to_csv(OUT / 'CTRP_treatment_identity_and_combination_flags.csv', index=False)
    total = finite = invalid_dose = below0 = above1 = 0
    single_rows = 0
    pairset, experiment_ids, compound_ids = set(), set(), set()
    ranges = []
    # Whole-file inspection; aggregate experiment x compound rather than every dose row.
    for chunk in pd.read_csv(base / 'v20.data.per_cpd_post_qc.txt', sep='\t', chunksize=500000):
        total += len(chunk)
        single_rows += int((~chunk.master_cpd_id.isin(combination_ids)).sum())
        values = chunk[['cpd_avg_pv', 'cpd_pred_pv', 'cpd_conc_umol']].to_numpy()
        finite += int(np.isfinite(values).all(axis=1).sum())
        invalid_dose += int((~np.isfinite(chunk.cpd_conc_umol) | chunk.cpd_conc_umol.le(0)).sum())
        below0 += int(chunk.cpd_avg_pv.lt(0).sum()); above1 += int(chunk.cpd_avg_pv.gt(1).sum())
        ranges.append((float(chunk.cpd_avg_pv.min()), float(chunk.cpd_avg_pv.max())))
        pairset.update(map(tuple, chunk[['experiment_id', 'master_cpd_id']].drop_duplicates().to_numpy()))
        experiment_ids.update(chunk.experiment_id.unique()); compound_ids.update(chunk.master_cpd_id.unique())
    # The source repeats experiment IDs with differing dates. Cell identity is invariant;
    # reduce only the ID-to-cell mapping, never arbitrarily select experiment metadata.
    assert ex.groupby('experiment_id').master_ccl_id.nunique().le(1).all()
    ex_to_cell = ex[['experiment_id', 'master_ccl_id']].drop_duplicates()
    pairs = pd.DataFrame(sorted(pairset), columns=['experiment_id', 'master_cpd_id']).merge(ex_to_cell, on='experiment_id', how='left', validate='many_to_one')
    RESULT['CTRPv2_2015'] = {'dose_rows': total, 'rows_finite_measured_fitted_dose': finite,
        'invalid_nonpositive_dose_rows': invalid_dose, 'experiment_compound_curves': len(pairs),
        'distinct_cell_drug_pairs': len(pairs.drop_duplicates(['master_ccl_id', 'master_cpd_id'])),
        'cells_in_measured_table': pairs.master_ccl_id.nunique(), 'compounds_in_measured_table': len(compound_ids),
        'cell_metadata_rows': len(cells), 'compound_metadata_rows': len(compounds),
        'experiment_metadata_rows': len(ex), 'unique_experiment_ids': ex.experiment_id.nunique(),
        'experiment_ids_with_repeated_metadata_rows': int(ex.experiment_id.value_counts().gt(1).sum()),
        'explicit_fixed_ratio_combination_entries': int(combo.sum()),
        'entries_without_explicit_combination_name': int((~combo).sum()),
        'dose_rows_without_explicit_combination_name': single_rows,
        'dose_rows_explicit_combination': total - single_rows,
        'distinct_cell_drug_pairs_without_explicit_combination': len(pairs[~pairs.master_cpd_id.isin(combination_ids)].drop_duplicates(['master_ccl_id','master_cpd_id'])),
        'measured_relative_value_range': [min(x[0] for x in ranges), max(x[1] for x in ranges)],
        'measured_values_below_0': below0, 'measured_values_above_1': above1,
        'qualification': 'cpd_avg_pv is aggregated measured relative viability on ratio scale; cpd_pred_pv is curve fit. Values above 1 need not be errors; do not divide by 100 or automatically clip.'}


def prism():
    p = DATA / 'raw/PRISM_19Q4/secondary-screen-replicate-collapsed-logfold-change.csv'
    d = pd.read_csv(p, index_col=0)
    a = d.to_numpy()
    RESULT['PRISM_secondary'] = {'matrix_shape': list(a.shape), 'finite_measurements': int(np.isfinite(a).sum()),
        'missing_measurements': int(np.isnan(a).sum()), 'finite_fraction': float(np.isfinite(a).mean()),
        'LFC_distribution': distribution(a), 'qualification': 'Measured barcode abundance after screening, relative to vehicle, with batch correction; independent of RNA assay, but not a single-cell death label. Fitted AUC/IC50 are separate.'}


def snapshot():
    status = json.loads((DATA / 'download_status.json').read_text())
    done = [f for f in status['files'] if f['status'] == '已完成']
    mismatch = []
    for f in done:
        p = DATA / f['path']
        if not p.exists() or p.stat().st_size != f['bytes']: mismatch.append(f['path'])
    RESULT['download_snapshot'] = {'updated_at': status['updated_at'], 'complete': len(done),
        'total': len(status['files']), 'complete_bytes': sum(f['bytes'] for f in done),
        'completed_file_size_or_presence_mismatches': mismatch,
        'pending': [{k: f.get(k) for k in ['dataset', 'name', 'status', 'downloaded_bytes', 'expected_size', 'error']} for f in status['files'] if f['status'] != '已完成']}


if __name__ == '__main__':
    snapshot(); save()
    print('HMS measurements', flush=True)
    cond, keys = functional(); save()
    print('Actual baseline and structures', flush=True)
    genes, baseline_cells, structures = baseline_and_structures(cond); save()
    print('L5 measured genes and L3 control coverage', flush=True)
    expression(cond, keys, genes, baseline_cells, structures); save()
    print('sciPlex actual conditions', flush=True)
    sciplex(); save()
    print('CTRP whole measured table', flush=True)
    ctrp(); save()
    print('PRISM measured matrix', flush=True)
    prism(); snapshot(); save()
    print('Finished: ' + str(OUT), flush=True)
