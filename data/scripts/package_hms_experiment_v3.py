#!/usr/bin/env python3
"""Reproducible experimental packaging of audited V2; no biological-QC release.

Run from any directory. Keeps raw/V1/V2 unchanged. No models are fitted here.
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'data/integration_20261106/v2_库板来源核对'
OUT = ROOT / 'data/prepared_20261106/HMS_v3_实验准备'
SEED = 20261106


def digest(s):
    return hashlib.sha256(str(s).encode()).hexdigest()


def balanced_folds(frame, group):
    # Assign whole identity groups, balancing distinct cell/drug/dose conditions.
    counts = frame.drop_duplicates('biological_condition_group').groupby(group).size()
    loads = [0] * 5
    mapping = {}
    for name, size in sorted(counts.items(), key=lambda x: (-x[1], digest(f'{SEED}|{x[0]}'))):
        fold = min(range(5), key=lambda k: (loads[k], k))
        mapping[name] = fold
        loads[fold] += int(size)
    return mapping


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(SRC / '核心库板条件表.csv').sort_values('library_condition_id').reset_index(drop=True)
    sig = pd.read_csv(SRC / '核心表达签名与矩阵行号.csv')
    z = np.load(SRC / '核心表达签名矩阵.npz', allow_pickle=False)
    base = pd.read_csv(SRC / 'CCLE背景矩阵_与表达基因同顺序.csv', index_col=0)
    bio = pd.read_csv(SRC / '库板内生物重复与技术孔引用.csv')
    assert c.library_condition_id.is_unique and sig.sig_id.is_unique
    assert np.array_equal(sig.sig_id.to_numpy(dtype=str), z['sig_id'])
    assert np.array_equal(base.columns.to_numpy(dtype=str), z['gene_id'])
    assert sig.groupby('library_condition_id').time_h.agg(set).apply(lambda x: x == {3., 24.}).all()
    assert len(sig) == 2 * len(c)

    chem = c[['pert_id', 'pert_iname', 'canonical_smiles']].drop_duplicates().sort_values('pert_id').reset_index(drop=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fingerprints = []
    for i, r in chem.iterrows():
        original = Chem.MolFromSmiles(r.canonical_smiles)
        assert original is not None
        parent = rdMolStandardize.Uncharger().uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(original)))
        stereo = Chem.MolToSmiles(parent, isomericSmiles=True)
        connectivity = Chem.MolToSmiles(parent, isomericSmiles=False)
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=parent, includeChirality=False)
        chem.loc[i, 'standardized_isomeric_smiles'] = stereo
        chem.loc[i, 'connectivity_smiles_for_split_only'] = connectivity
        chem.loc[i, 'chemical_split_group'] = 'CHEM_' + digest(connectivity)[:16]
        chem.loc[i, 'scaffold_smiles'] = scaffold
        chem.loc[i, 'scaffold_split_group'] = 'SCAF_' + digest(scaffold or 'ACYCLIC')[:16]
        chem.loc[i, 'original_fragment_count'] = len(Chem.GetMolFrags(original))
        fingerprints.append(generator.GetFingerprintAsNumPy(parent))
    assert chem.pert_id.is_unique
    chem['drug_index'] = np.arange(len(chem))
    c = c.merge(chem[['pert_id', 'drug_index', 'chemical_split_group', 'scaffold_split_group']], on='pert_id', validate='many_to_one')
    cells = c[['cell_id', 'ModelID']].drop_duplicates().sort_values('ModelID').reset_index(drop=True)
    cells['cell_index'] = np.arange(len(cells))
    c = c.merge(cells, on=['cell_id', 'ModelID'], validate='many_to_one')
    c['sample_index'] = np.arange(len(c))
    c['n_library_views_of_condition'] = c.groupby('biological_condition_group').library_condition_id.transform('size')
    c['condition_weight'] = 1.0 / c.n_library_views_of_condition
    c['eligibility_stage'] = 'numeric_source_checks_passed; biological_QC_pending'
    c['cell_drug_split_group'] = c.ModelID + '|' + c.chemical_split_group
    assert np.allclose(c.groupby('biological_condition_group').condition_weight.sum(), 1)

    lookup = dict(zip(z['sig_id'].astype(str), range(len(z['sig_id']))))
    expression = np.empty((len(c), 2, len(z['gene_id'])), dtype=np.float32)
    signature_ids = np.empty((len(c), 2), dtype=z['sig_id'].dtype)
    for j, t in enumerate([3, 24]):
        ids = c[f'signature_ids_{t}h'].astype(str)
        assert not ids.str.contains('|', regex=False).any(), 'Multiple signatures require explicit aggregation.'
        expression[:, j] = z['expression_z'][[lookup[s] for s in ids]]
        signature_ids[:, j] = ids
    selected = bio[bio.library_condition_id.isin(c.library_condition_id)]
    bio_files = sorted(selected.source_file.unique())
    bv = selected.pivot(index='library_condition_id', columns='source_file', values='relative_cell_count').reindex(index=c.library_condition_id, columns=bio_files).to_numpy()
    assert bv.shape == (len(c), 3) and np.isfinite(bv).all()
    assert np.allclose(bv.mean(1), c.relative_cell_count)
    assert np.allclose(bv.std(1, ddof=1), c.relative_cell_count_sd)
    background = base.loc[cells.ModelID].to_numpy(dtype=np.float32)
    assert np.isfinite(expression).all() and np.isfinite(background).all()
    np.savez_compressed(OUT / 'arrays.npz',
        expression_z=expression, expression_time_h=np.array([3, 24], dtype=np.float32),
        expression_signature_id=signature_ids, expression_mask=np.ones((len(c), 2), dtype=bool),
        relative_cell_count_72h=c.relative_cell_count.to_numpy(),
        relative_net_growth_72h=c.relative_net_growth.to_numpy(),
        count_bioreplicates=bv, bioreplicate_source=np.asarray(bio_files, dtype=str),
        functional_mask=np.ones(len(c), dtype=bool),
        baseline_log2TPMp1=background, gene_id=z['gene_id'].astype(str),
        drug_morgan=np.stack(fingerprints).astype(np.uint8),
        cell_index=c.cell_index.to_numpy(), drug_index=c.drug_index.to_numpy(),
        log10_dose_uM=np.log10(c.dose_uM.to_numpy()).astype(np.float32),
        condition_weight=c.condition_weight.to_numpy(),
        library_condition_id=c.library_condition_id.to_numpy(dtype=str))
    c.to_csv(OUT / 'samples.csv', index=False)
    chem.to_csv(OUT / 'drugs.csv', index=False)
    cells.to_csv(OUT / 'cells.csv', index=False)
    selected.to_csv(OUT / 'functional_replicates.csv', index=False)
    sig.to_csv(OUT / 'expression_provenance.csv', index=False)
    genes = pd.read_csv(ROOT / 'data/raw/LINCS_GSE70138/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz', sep='\t').set_index('pr_gene_id')
    genes.loc[z['gene_id'].astype(int)].reset_index().to_csv(OUT / 'genes.csv', index=False)

    # Draft split manifests. All times/libraries/repeats of a held-out unit stay together.
    splits = []
    for protocol, group in [('unseen_drug', 'chemical_split_group'), ('unseen_scaffold', 'scaffold_split_group'), ('unseen_cell', 'ModelID'), ('unseen_cell_drug_pair', 'cell_drug_split_group')]:
        if protocol == 'unseen_cell_drug_pair':
            dr = {g: i for i, g in enumerate(sorted(c.chemical_split_group.unique()))}
            cl = {g: i for i, g in enumerate(sorted(c.ModelID.unique()))}
            assignment = (c.chemical_split_group.map(dr) + c.ModelID.map(cl)) % 5
        else:
            assignment = c[group].map(balanced_folds(c, group))
        for fold in range(5):
            role = np.where(assignment.eq(fold), 'test', np.where(assignment.eq((fold + 1) % 5), 'validation', 'train'))
            d = c[['sample_index', 'library_condition_id', 'biological_condition_group', group]].copy().rename(columns={group: 'heldout_group'})
            d['protocol'] = protocol
            d['fold'] = fold
            d['role'] = role
            assert d.groupby('heldout_group').role.nunique().eq(1).all()
            if protocol == 'unseen_cell_drug_pair':
                for column in ['ModelID', 'chemical_split_group']:
                    assert set(c.loc[role != 'train', column]) <= set(c.loc[role == 'train', column])
            splits.append(d)
    for protocol, test_dose, val_dose in [('dose_interpolation', .37, 1.11), ('dose_extrapolation', 10., 3.33)]:
        d = c[['sample_index', 'library_condition_id', 'biological_condition_group']].copy()
        d['heldout_group'] = c.biological_condition_group
        d['protocol'] = protocol
        d['fold'] = 0
        d['role'] = np.where(np.isclose(c.dose_uM, test_dose), 'test', np.where(np.isclose(c.dose_uM, val_dose), 'validation', 'train'))
        train = c.loc[d.role.eq('train')].groupby('cell_drug_split_group').dose_uM.agg(['min', 'max'])
        for ix in d.index[d.role.ne('train')]:
            r = c.loc[ix]
            limits = train.loc[r.cell_drug_split_group]
            if protocol == 'dose_interpolation':
                assert limits['min'] < r.dose_uM < limits['max']
            else:
                assert limits['max'] < r.dose_uM
        splits.append(d)
    splits = pd.concat(splits, ignore_index=True)
    assert splits.groupby(['protocol', 'fold', 'biological_condition_group']).role.nunique().eq(1).all()
    splits.to_csv(OUT / 'split_assignments_DRAFT.csv', index=False)
    summary = splits.groupby(['protocol', 'fold', 'role']).agg(library_records=('library_condition_id', 'size'), distinct_conditions=('biological_condition_group', 'nunique')).reset_index()
    summary.to_csv(OUT / 'split_counts.csv', index=False)

    comparisons = []
    for group, rows in c.groupby('biological_condition_group'):
        if len(rows) == 2:
            a, b = rows.index
            r = {'biological_condition_group': group, 'cell_id': c.loc[a, 'cell_id'], 'pert_id': c.loc[a, 'pert_id'],
                 'pert_iname': c.loc[a, 'pert_iname'], 'dose_uM': c.loc[a, 'dose_uM'],
                 'library_a': c.loc[a, 'library_id'], 'library_b': c.loc[b, 'library_id'],
                 'absolute_count_difference': abs(c.loc[a, 'relative_cell_count'] - c.loc[b, 'relative_cell_count'])}
            for j, t in enumerate([3, 24]):
                r[f'expression_pearson_{t}h'] = float(np.corrcoef(expression[a, j], expression[b, j])[0, 1])
            comparisons.append(r)
    comparisons = pd.DataFrame(comparisons)
    comparisons.to_csv(OUT / 'cross_library_descriptive_checks.csv', index=False)
    def quantiles(x):
        return {str(q): float(np.quantile(x, q)) for q in [0, .25, .5, .75, .95, 1]}
    p1 = pd.read_csv(ROOT / 'data/raw/LINCS_GSE92742/GSE92742_Broad_LINCS_gene_info.txt.gz', sep='\t')
    delta = pd.read_csv(ROOT / 'data/raw/LINCS_GSE92742/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt.gz', sep='\t')
    lm1 = set(p1.loc[p1.pr_is_lm.eq(1), 'pr_gene_id'])
    lm2 = set(genes.index[genes.pr_is_lm.eq(1)])
    assert lm1 == lm2 == set(z['gene_id'].astype(int))
    metrics = sig[['distil_nsample', 'tas', 'distil_cc_q75', 'distil_ss']].apply(pd.to_numeric, errors='coerce')
    listed = sig.distil_id.str.split('|').str.len()
    assert np.array_equal(listed.to_numpy(), metrics.distil_nsample.to_numpy())
    assert np.array_equal(metrics.distil_cc_q75.isna().to_numpy(), metrics.distil_nsample.eq(1).to_numpy())
    audit = {
        'stage': 'experiment_ready_arrays_and_DRAFT_splits; biological_QC_pending',
        'models_trained': 0, 'final_training_release': False,
        'library_records': len(c), 'distinct_conditions': c.biological_condition_group.nunique(),
        'cell_drug_pairs': c.cell_drug_split_group.nunique(), 'cells': len(cells), 'drug_ids': len(chem),
        'chemical_connectivity_groups': chem.chemical_split_group.nunique(), 'scaffold_groups': chem.scaffold_split_group.nunique(),
        'expression_shape': list(expression.shape), 'baseline_shape': list(background.shape),
        'sum_condition_weights': float(c.condition_weight.sum()),
        'conditions_with_two_libraries': len(comparisons),
        'functional_bioreplicate_rows': len(selected),
        'technical_wells_in_core': int(selected.n_technical_rows.sum()),
        'technical_rows_per_bioreplicate': selected.n_technical_rows.value_counts().to_dict(),
        'functional_bioreplicate_spearman_descriptive': {f'{i+1}_{j+1}': float(spearmanr(bv[:, i], bv[:, j]).statistic) for i in range(3) for j in range(i+1, 3)},
        'count_sd_across_bioreplicates_quantiles': quantiles(c.relative_cell_count_sd),
        'cross_library_count_difference_quantiles': quantiles(comparisons.absolute_count_difference),
        'cross_library_expression_pearson_quantiles': {str(t): quantiles(comparisons[f'expression_pearson_{t}h']) for t in [3, 24]},
        'expression_metric_missing_counts': metrics.isna().sum().to_dict(),
        'missing_replicate_correlation_all_and_only_single_instance_signatures': True,
        'expression_replicate_number_distribution': metrics.distil_nsample.value_counts().sort_index().to_dict(),
        'phaseI_phaseII_pr_is_lm_genes_equal': True,
        'delta_landmark_file_intersection_with_measured_landmarks': len(set(delta.pr_gene_id) & lm2),
        'split_assertions_passed': ['all_library_views_and_times_grouped', 'heldout_identity_disjoint', 'pair_test_entities_present_in_train', 'interpolation_within_training_dose_range', 'extrapolation_above_training_dose_range'],
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'rdkit': rdkit.__version__},
        'seed': SEED,
        'notes': ['Descriptive repeat statistics are not independent biological validation.', 'Connectivity grouping is for conservative splits only; original chemical identities and stereochemistry are retained.', 'No TAS/response-size-based exclusion and no clipping of function labels.']}
    def default(x):
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    (OUT / 'audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2, default=default) + '\n')
    files = []
    for path in sorted(SRC.glob('*')):
        if path.is_file():
            files.append({'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    (OUT / 'input_fingerprints.json').write_text(json.dumps(files, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(audit, ensure_ascii=False, indent=2, default=default))


if __name__ == '__main__':
    main()
