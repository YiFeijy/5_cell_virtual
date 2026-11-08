#!/usr/bin/env python3
"""Outer-test metrics, grouped errors, ranking and paired drug bootstrap."""
import os
for v in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[v] = '2'
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from train import OUT, load, METHODS, SEEDS


def avg(values, weights):
    valid = np.isfinite(values)
    return float(np.average(values[valid], weights=weights[valid])) if valid.any() else np.nan


def row_errors(indices, pe, pf):
    a, _, _ = load()
    errors = {}
    if pe is not None:
        truth = a['expression_z24'][indices].astype(np.float64)
        delta = pe-truth
        errors['expression_mse'] = (delta**2).mean(1)
        errors['expression_mae'] = abs(delta).mean(1)
        x, y = pe-pe.mean(1, keepdims=True), truth-truth.mean(1, keepdims=True)
        den = np.sqrt((x*x).sum(1)*(y*y).sum(1))
        errors['expression_profile_pearson'] = np.divide((x*y).sum(1), den, out=np.full(len(indices), np.nan), where=den>1e-12)
    if pf is not None:
        delta = pf-a['count72'][indices]
        errors['function_mae'] = abs(delta)
        errors['function_mse'] = delta**2
    return errors


def metrics(indices, pe, pf):
    a, s, _ = load()
    frame = s.iloc[indices]
    weights = 1./frame.groupby('biological_condition_group').library_condition_id.transform('size').to_numpy()
    errors = row_errors(indices, pe, pf)
    result = {'library_records': len(indices), 'conditions': frame.biological_condition_group.nunique(),
              'drugs': frame.pert_id.nunique(), 'cells': frame.cell_id.nunique()}
    result.update({k: avg(v, weights) for k, v in errors.items()})
    if pf is not None:
        result['function_rmse'] = result['function_mse']**.5
        yf = a['count72'][indices]
        variance = np.average((yf-np.average(yf, weights=weights))**2, weights=weights)
        result['function_r2'] = 1-result['function_mse']/variance if variance else np.nan
        group = pd.DataFrame({'group': frame.biological_condition_group.to_numpy(), 'truth': yf, 'pred': pf}).groupby('group')[['truth', 'pred']].mean()
        result['function_condition_spearman'] = float(spearmanr(group.truth, group.pred).statistic) if group.pred.nunique()>1 else np.nan
        drug_error = pd.DataFrame({'drug': frame.pert_id.to_numpy(), 'err': errors['function_mae']*weights, 'weight': weights}).groupby('drug')[['err', 'weight']].sum()
        result['function_drug_macro_mae'] = float((drug_error.err/drug_error.weight).mean())
    if pe is not None:
        result['expression_valid_profile_correlations'] = int(np.isfinite(errors['expression_profile_pearson']).sum())
    return result


def combine(protocol, method, seed, baseline=False):
    folds = range(5) if protocol == 'unseen_drug' else [0]
    pieces = []
    for fold in folds:
        path = OUT / ('baselines' if baseline else 'runs') / (f'{protocol}__fold{fold}__{method}' if baseline else f'{protocol}__fold{fold}__{method}__seed{seed}')
        with np.load(path / 'test_predictions.npz', allow_pickle=False) as z:
            pieces.append({k: z[k] for k in z.files if k!='latent'})
    ids = np.concatenate([p['sample_index'] for p in pieces])
    assert len(ids)==len(set(ids))
    order = np.argsort(ids)
    pe = np.concatenate([p['expression_predicted'] for p in pieces])[order] if 'expression_predicted' in pieces[0] else None
    pf = np.concatenate([p['count72_predicted'] for p in pieces])[order] if 'count72_predicted' in pieces[0] else None
    return ids[order], pe, pf


def ranking(protocol, method, seed, ix, pf):
    if pf is None:
        return []
    a, s, sp = load()
    foldmap = sp[sp.protocol.eq(protocol) & sp.role.eq('test')].set_index('sample_index').fold
    frame = s.iloc[ix][['biological_condition_group', 'cell_id', 'pert_id', 'dose_uM']].copy()
    frame['fold'] = foldmap.reindex(ix).to_numpy()
    frame['predicted'] = pf
    frame['observed'] = a['count72'][ix]
    # Rank drugs measured by the SAME fitted fold model; do not mix fold calibrations.
    frame = frame.groupby(['fold', 'cell_id', 'pert_id', 'dose_uM'], as_index=False)[['predicted', 'observed']].mean()
    rows = []
    for (fold, cell, dose), group in frame.groupby(['fold', 'cell_id', 'dose_uM']):
        if len(group)<10:
            continue
        # Equal predictions have no preferred drug: exact expected random tie break.
        oracle = group.nsmallest(5, 'observed')
        ties = group.predicted.nunique()<len(group)
        cutoff = np.sort(group.predicted.to_numpy())[4]
        must = group[group.predicted<cutoff]
        tied = group[group.predicted.eq(cutoff)]
        remaining = 5-len(must)
        chosen = float((must.observed.sum()+remaining*tied.observed.mean())/5)
        overlap = (len(set(must.pert_id)&set(oracle.pert_id))+remaining*len(set(tied.pert_id)&set(oracle.pert_id))/len(tied))/5
        rows.append({'protocol': protocol, 'method': method, 'seed': seed, 'fold': fold, 'cell_id': cell, 'dose_uM': dose,
                     'n_drugs': len(group), 'top5_observed_count': chosen, 'oracle_top5_count': float(oracle.observed.mean()),
                     'all_drugs_mean_count': float(group.observed.mean()), 'top5_regret': chosen-float(oracle.observed.mean()),
                     'top5_overlap': float(overlap), 'ties_averaged_exactly': ties})
    return rows


def main():
    reports = OUT / 'metrics'; reports.mkdir(exist_ok=True)
    dest = OUT / 'predictions'; dest.mkdir(exist_ok=True)
    a, s, sp = load()
    progress = pd.read_csv(OUT / 'training_runs.csv')
    assert len(progress)==144
    baseline_names = ['zero_response', 'train_mean', 'cell_mean', 'cell_dose_mean', 'ridge_independent']
    all_metrics, all_grouped, all_ranking, all_drug = [], [], [], []
    sample_errors = {}
    for protocol in ['unseen_drug', 'dose_extrapolation']:
        methods = baseline_names+METHODS
        for method in methods:
            baseline = method in baseline_names
            seeds = ['deterministic'] if baseline else SEEDS
            runs = []
            for seed in seeds:
                ix, pe, pf = combine(protocol, method, seed, baseline)
                runs.append((ix, pe, pf))
            if not baseline:
                ix = runs[0][0]
                assert all(np.array_equal(ix, r[0]) for r in runs)
                runs.append((ix, np.mean([r[1] for r in runs], axis=0) if runs[0][1] is not None else None,
                             np.mean([r[2] for r in runs], axis=0) if runs[0][2] is not None else None))
                seeds = list(seeds)+['ensemble']
            for seed, (ix, pe, pf) in zip(seeds, runs):
                payload = {'sample_index': ix}
                if pe is not None: payload['expression_predicted'] = pe
                if pf is not None: payload['count72_predicted'] = pf
                np.savez_compressed(dest / f'{protocol}__{method}__{seed}.npz', **payload)
                subsets = {'all': np.ones(len(ix), bool), 'rep3plus': s.iloc[ix].sensitivity_rep3plus.to_numpy(),
                           'without_CGP60474': s.iloc[ix].sensitivity_without_CGP60474.to_numpy()}
                for scope, mask in subsets.items():
                    all_metrics.append({'protocol': protocol, 'method': method, 'seed': seed, 'scope': scope,
                                        **metrics(ix[mask], pe[mask] if pe is not None else None, pf[mask] if pf is not None else None)})
                all_ranking.extend(ranking(protocol, method, seed, ix, pf))
                if seed in ['ensemble', 'deterministic']:
                    frame = s.iloc[ix]
                    for column in ['cell_id', 'dose_uM']:
                        for group, rows in frame.groupby(column):
                            positions = np.flatnonzero(frame.index.isin(rows.index))
                            all_grouped.append({'protocol': protocol, 'method': method, 'group_type': column, 'group': group,
                                                **metrics(ix[positions], pe[positions] if pe is not None else None, pf[positions] if pf is not None else None)})
            # Inference uncertainty uses average SINGLE-SEED errors, not ensemble errors.
            errors = [row_errors(r[0], r[1], r[2]) for r in runs[:1 if baseline else 3]]
            means = {key: np.mean([e[key] for e in errors], axis=0) for key in errors[0]}
            sample_errors[(protocol, method)] = (runs[0][0], means)
            ix = runs[0][0]
            frame = s.iloc[ix]
            for drug, rows in frame.groupby('pert_id'):
                positions = np.flatnonzero(frame.index.isin(rows.index))
                w = a['condition_weight'][ix[positions]]
                all_drug.append({'protocol': protocol, 'method': method, 'pert_id': drug, 'drug_name': rows.pert_iname.iloc[0],
                                 'conditions': rows.biological_condition_group.nunique(), **{k: avg(v[positions], w) for k,v in means.items()}})
    detail = pd.DataFrame(all_metrics)
    detail.to_csv(reports / 'all_metrics.csv', index=False)
    selected = detail[detail.seed.astype(str).ne('ensemble')]
    metric_names = [c for c in detail if c not in ['protocol', 'method', 'seed', 'scope', 'library_records', 'conditions', 'drugs', 'cells', 'expression_valid_profile_correlations']]
    summary = selected.groupby(['protocol', 'method', 'scope'])[metric_names].agg(['mean', 'std', 'count'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary.reset_index().to_csv(reports / 'seed_summary.csv', index=False)
    pd.DataFrame(all_grouped).to_csv(reports / 'grouped_metrics.csv', index=False)
    pd.DataFrame(all_drug).to_csv(reports / 'per_drug_errors.csv', index=False)
    ranks = pd.DataFrame(all_ranking)
    ranks.to_csv(reports / 'ranking_by_fold_cell_dose.csv', index=False)
    ranks.groupby(['protocol', 'method', 'seed'])[['top5_observed_count','top5_regret','top5_overlap']].mean().reset_index().to_csv(reports / 'ranking_seed_summary.csv', index=False)
    differences = []
    for protocol in ['unseen_drug', 'dose_extrapolation']:
        for method in ['joint_shared', 'joint_shared_E01', 'joint_shared_F01', 'joint_private', 'joint_bottleneck']:
            for metric, reference in [('function_mae', 'func_only'), ('expression_mse', 'expr_only')]:
                ix, values = sample_errors[(protocol, method)]
                ref_ix, ref_values = sample_errors[(protocol, reference)]
                assert np.array_equal(ix, ref_ix)
                delta = values[metric]-ref_values[metric]
                w = a['condition_weight'][ix]
                group = pd.DataFrame({'drug': s.iloc[ix].pert_id.to_numpy(), 'numerator': delta*w, 'mass': w}).groupby('drug')[['numerator','mass']].sum()
                rng = np.random.default_rng(20261110)
                bootstrap = rng.integers(0,len(group),(2000,len(group)))
                b = group.numerator.to_numpy()[bootstrap].sum(1)/group.mass.to_numpy()[bootstrap].sum(1)
                differences.append({'protocol':protocol,'method':method,'reference':reference,'metric':metric,
                                    'difference_mean':float(np.average(delta,weights=w)), 'ci95_low':float(np.quantile(b,.025)),
                                    'ci95_high':float(np.quantile(b,.975)), 'resampling_unit':'drug', 'n_drugs':len(group),
                                    'seed_treatment':'mean per-observation error across 3 seeds; no ensemble averaging of predictions'})
        ix, values = sample_errors[(protocol,'joint_bottleneck')]
        _, refs = sample_errors[(protocol,'function_bottleneck_control')]
        delta = values['function_mae']-refs['function_mae']; w=a['condition_weight'][ix]
        group=pd.DataFrame({'drug':s.iloc[ix].pert_id.to_numpy(),'numerator':delta*w,'mass':w}).groupby('drug')[['numerator','mass']].sum()
        rng=np.random.default_rng(20261110);bix=rng.integers(0,len(group),(2000,len(group)))
        b=group.numerator.to_numpy()[bix].sum(1)/group.mass.to_numpy()[bix].sum(1)
        differences.append({'protocol':protocol,'method':'joint_bottleneck','reference':'function_bottleneck_control','metric':'function_mae',
                            'difference_mean':float(np.average(delta,weights=w)),'ci95_low':float(np.quantile(b,.025)),'ci95_high':float(np.quantile(b,.975)),
                            'resampling_unit':'drug','n_drugs':len(group),'seed_treatment':'mean per-observation error across 3 seeds'})
    pd.DataFrame(differences).to_csv(reports / 'paired_drug_bootstrap.csv',index=False)
    print(summary.reset_index().query("protocol=='unseen_drug' and scope=='all'")[[('method'),('function_mae_mean'),('expression_mse_mean')]].to_string(index=False),flush=True)
    print(pd.DataFrame(differences).query("protocol=='unseen_drug'").to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
