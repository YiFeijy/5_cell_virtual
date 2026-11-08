#!/usr/bin/env python3
"""Independent Ridge baselines for the HMS preparation bundle.

Default validates loading/split/identity only. --fit fits ONE task and ONE fold,
selecting alpha on validation, then exports test predictions. This is a pilot
on pre-QC data, never a final study result. No expression label enters inputs.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/prepared_20261106/HMS_v3_实验准备'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fit', action='store_true')
    p.add_argument('--protocol', default='unseen_drug', choices=['unseen_drug', 'unseen_scaffold', 'unseen_cell', 'unseen_cell_drug_pair', 'dose_interpolation', 'dose_extrapolation'])
    p.add_argument('--fold', type=int, default=0)
    p.add_argument('--task', default='function', choices=['function', 'expression3', 'expression24'])
    a = p.parse_args()
    z = np.load(DATA / 'arrays.npz', allow_pickle=False)
    s = pd.read_csv(DATA / 'samples.csv')
    d = pd.read_csv(DATA / 'split_assignments_DRAFT.csv')
    d = d[d.protocol.eq(a.protocol) & d.fold.eq(a.fold)].sort_values('sample_index')
    assert len(d) == len(s), 'Protocol/fold not present.'
    assert np.array_equal(d.library_condition_id.to_numpy(dtype=str), z['library_condition_id'])
    assert np.array_equal(s.library_condition_id.to_numpy(dtype=str), z['library_condition_id'])
    x = np.concatenate([z['baseline_log2TPMp1'][z['cell_index']], z['drug_morgan'][z['drug_index']], z['log10_dose_uM'][:, None]], axis=1).astype(np.float64)
    y = z['relative_cell_count_72h'] if a.task == 'function' else z['expression_z'][:, 0 if a.task == 'expression3' else 1]
    roles = {r: d.role.eq(r).to_numpy() for r in ['train', 'validation', 'test']}
    assert np.isfinite(x).all() and np.isfinite(y).all()
    assert d.groupby('biological_condition_group').role.nunique().eq(1).all()
    for r in roles:
        assert roles[r].any()
    report = {'stage': 'pre_QC_pilot', 'protocol': a.protocol, 'fold': a.fold, 'task': a.task,
              'input_shape': list(x.shape), 'target_shape': list(y.shape),
              'role_library_records': {r: int(m.sum()) for r, m in roles.items()},
              'checks': ['array_ID_order', 'finite_inputs_targets', 'no_condition_group_crosses_roles', 'no_observed_expression_as_input']}
    if not a.fit:
        report['models_fitted_this_run'] = 0
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return
    tr, va, te = (roles[r] for r in ['train', 'validation', 'test'])
    w = z['condition_weight']
    scaler = StandardScaler().fit(x[tr], sample_weight=w[tr])
    xt = scaler.transform(x)
    trials = []
    best = None
    for alpha in [1., 10., 100., 1000.]:
        model = Ridge(alpha=alpha, solver='lsqr').fit(xt[tr], y[tr], sample_weight=w[tr])
        pred = model.predict(xt[va])
        err = (pred - y[va]) ** 2
        mse = float(np.average(err if err.ndim == 1 else err.mean(1), weights=w[va]))
        trials.append({'alpha': alpha, 'validation_mse': mse})
        if best is None or mse < best[0]:
            best = (mse, model, alpha)
    prediction = best[1].predict(xt[te])
    err = prediction - y[te]
    per_mse = err**2 if err.ndim == 1 else (err**2).mean(1)
    per_mae = np.abs(err) if err.ndim == 1 else np.abs(err).mean(1)
    rows = s.loc[te, ['library_condition_id', 'biological_condition_group', 'cell_id', 'pert_id', 'dose_uM', 'condition_weight']].copy()
    rows['mse'] = per_mse
    rows['mae'] = per_mae
    if prediction.ndim == 1:
        rows['observed'] = y[te]
        rows['predicted'] = prediction
        report['test_spearman_descriptive'] = float(spearmanr(y[te], prediction).statistic)
    else:
        rows['profile_pearson'] = [float(np.corrcoef(u, v)[0, 1]) if np.std(u) > 0 and np.std(v) > 0 else np.nan for u, v in zip(y[te], prediction)]
    report.update({'models_fitted_this_run': len(trials), 'alpha_trials': trials, 'chosen_alpha': best[2],
                   'test_mse': float(np.average(per_mse, weights=w[te])),
                   'test_mae': float(np.average(per_mae, weights=w[te])),
                   'interpretation': 'One pilot fold; final biological eligibility and confirmatory experiment protocol remain separate.'})
    out = ROOT / 'data/baseline_pilots_20261106' / f'{a.protocol}_fold{a.fold}_{a.task}'
    out.mkdir(parents=True, exist_ok=True)
    rows.to_csv(out / 'test_observation_errors.csv', index=False)
    np.savez_compressed(out / 'test_predictions.npz', predicted=prediction, observed=y[te], library_condition_id=z['library_condition_id'][te])
    (out / 'metrics.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
