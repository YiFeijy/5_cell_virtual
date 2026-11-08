#!/usr/bin/env python3
"""Training-only mean baselines and exact row-space weighted Ridge."""
import os
for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[v] = '2'
from pathlib import Path
import json
import numpy as np
import pandas as pd
from train import OUT, inputs, load


def main():
    a, s, _ = load()
    traces, checks = [], []
    for protocol, fold in [('unseen_drug', f) for f in range(5)] + [('dose_extrapolation', 0)]:
        xx, _, _, w, roles, _ = inputs(protocol, fold)
        tr, va = roles['train'], roles['validation']
        ye, yf = a['expression_z24'].astype(np.float64), a['count72']
        em, fm = np.average(ye[tr], axis=0, weights=w[tr]), float(np.average(yf[tr], weights=w[tr]))
        means_cell, means_dose = {}, {}
        for key, group in s.iloc[tr].groupby('cell_id'):
            ix = group.index.to_numpy()
            means_cell[key] = (np.average(ye[ix], axis=0, weights=w[ix]), np.average(yf[ix], weights=w[ix]))
        for key, group in s.iloc[tr].groupby(['cell_id', 'dose_uM']):
            ix = group.index.to_numpy()
            means_dose[key] = (np.average(ye[ix], axis=0, weights=w[ix]), np.average(yf[ix], weights=w[ix]))
        X = np.concatenate(xx, axis=1).astype(np.float64)
        xm = np.average(X[tr], axis=0, weights=w[tr])
        xc = X-xm
        anchors = set()
        for col in ['pert_id', 'cell_id', 'dose_uM']:
            anchors.update(s.iloc[tr].groupby(col).head(1).index)
        anchors = sorted(anchors)
        for attempt in range(120):
            _, singular, vt = np.linalg.svd(xc[anchors], full_matrices=False)
            keep = singular > max(singular[0]*1e-10, 1e-10)
            basis = vt[keep]
            reduced = xc @ basis.T
            residual = np.max(abs(xc[tr]-reduced[tr]@basis), axis=1)
            if residual.max() < 1e-7:
                break
            anchors.append(int(tr[np.argmax(residual)]))
        assert residual.max() < 1e-7
        # Orthonormal reduction preserves the train-row span and Ridge penalty.
        R = reduced
        gram = R[tr].T @ (w[tr, None]*R[tr])
        rhs_e = R[tr].T @ (w[tr, None]*(ye[tr]-em))
        rhs_f = R[tr].T @ (w[tr]*(yf[tr]-fm))
        choices = {}
        for task, rhs, truth, center in [('expression', rhs_e, ye, em), ('function', rhs_f, yf, fm)]:
            best = None
            for alpha in [1., 10., 100., 1000.]:
                coef = np.linalg.solve(gram+alpha*np.eye(len(gram)), rhs)
                prediction = R[va] @ coef+center
                err = (prediction-truth[va])**2
                score = float(np.average(err.mean(1) if err.ndim == 2 else err, weights=w[va]))
                traces.append({'protocol': protocol, 'fold': fold, 'task': task, 'alpha': alpha, 'validation_mse': score})
                if best is None or score < best[0]:
                    best = score, alpha, coef
            choices[task] = best
        check = {'protocol': protocol, 'fold': fold, 'training_feature_rank': len(basis), 'train_reconstruction_max_error': float(residual.max()),
                 'ridge_expression_alpha': choices['expression'][1], 'ridge_function_alpha': choices['function'][1]}
        if protocol == 'unseen_drug' and fold == 0:
            from sklearn.linear_model import Ridge
            # Independent numerical verification against the full feature solve.
            reference = Ridge(alpha=choices['function'][1], solver='lsqr', tol=1e-11).fit(X[tr], yf[tr], sample_weight=w[tr])
            exact = R[va] @ choices['function'][2]+fm
            difference = float(np.max(abs(exact-reference.predict(X[va]))))
            assert difference < 1e-6
            check['full_sklearn_ridge_prediction_max_difference'] = difference
        checks.append(check)
        for role, ix in [('validation', va), ('test', roles['test'])]:
            global_e = np.broadcast_to(em, (len(ix), len(em)))
            cell_values = [means_cell.get(r.cell_id, (em, fm)) for r in s.iloc[ix].itertuples()]
            dose_values = [means_dose.get((r.cell_id, r.dose_uM), means_cell.get(r.cell_id, (em, fm))) for r in s.iloc[ix].itertuples()]
            predictions = {
                'zero_response': (np.zeros_like(global_e), None),
                'train_mean': (global_e, np.full(len(ix), fm)),
                'cell_mean': (np.stack([v[0] for v in cell_values]), np.array([v[1] for v in cell_values])),
                'cell_dose_mean': (np.stack([v[0] for v in dose_values]), np.array([v[1] for v in dose_values])),
                'ridge_independent': (R[ix] @ choices['expression'][2]+em, R[ix] @ choices['function'][2]+fm)}
            for method, (pe, pf) in predictions.items():
                folder = OUT / 'baselines' / f'{protocol}__fold{fold}__{method}'
                folder.mkdir(parents=True, exist_ok=True)
                payload = {'sample_index': ix, 'expression_predicted': pe.astype(np.float32)}
                if pf is not None: payload['count72_predicted'] = pf
                np.savez_compressed(folder / f'{role}_predictions.npz', **payload)
                (folder / 'config.json').write_text(json.dumps(check, indent=2)+'\n')
        print(json.dumps(check), flush=True)
    pd.DataFrame(traces).to_csv(OUT / 'ridge_validation_search.csv', index=False)
    (OUT / 'baseline_numeric_checks.json').write_text(json.dumps(checks, indent=2)+'\n')


if __name__ == '__main__':
    main()
