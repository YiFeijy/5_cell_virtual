#!/usr/bin/env python3
"""Independent checks on split integrity, metrics, source values and artifacts."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error
from train import OUT,ROOT,load


def main():
    a,s,sp=load();frozen=json.loads((OUT/'frozen_training_protocol.json').read_text())
    assert hashlib.sha256((OUT/'data/arrays.npz').read_bytes()).hexdigest()==frozen['data_sha256']
    assert hashlib.sha256((OUT/'scripts/train.py').read_bytes()).hexdigest()==frozen['script_sha256']
    runs=pd.read_csv(OUT/'training_runs.csv');assert len(runs)==144
    expected={tuple(j) for j in frozen['jobs']}
    actual={(r.method,r.protocol,r.fold,r.seed) for r in runs.itertuples()};assert actual==expected
    assert not runs.test_used_for_checkpoint_selection.any()
    assert np.isfinite(a['expression_z24']).all() and np.isfinite(a['count72']).all()
    for (protocol,fold),d in sp.groupby(['protocol','fold']):
        assert d.groupby('biological_condition_group').role.nunique().eq(1).all()
        if protocol in ['unseen_drug','unseen_scaffold']:
            assert d.groupby('heldout_group').role.nunique().eq(1).all()
    dose=sp[sp.protocol.eq('dose_extrapolation')]
    tr=s.loc[dose[dose.role.eq('train')].sample_index].groupby(['cell_id','pert_id']).dose_uM.max()
    for r in s.loc[dose[dose.role.eq('test')].sample_index].itertuples():
        assert (r.cell_id,r.pert_id) in tr and tr.loc[(r.cell_id,r.pert_id)]<r.dose_uM
    checkpoint_checks=0
    for r in runs.itertuples():
        p=OUT/'runs'/f'{r.protocol}__fold{r.fold}__{r.method}__seed{r.seed}'
        h=pd.read_csv(p/'history.csv')
        assert h.loc[h.epoch.eq(r.selected_epoch),'validation_score'].iloc[0]<=h.validation_score.min()+1.1e-5
        assert (p/'model_state.pt').is_file()
        for role in ['validation','test']:
            z=np.load(p/f'{role}_predictions.npz',allow_pickle=False)
            ids=z['sample_index'];d=sp[sp.protocol.eq(r.protocol)&sp.fold.eq(r.fold)&sp.role.eq(role)]
            assert set(ids)==set(d.sample_index)
            for key in z.files:
                assert np.isfinite(z[key]).all()
        checkpoint_checks+=1
    z=np.load(OUT/'predictions/unseen_drug__joint_shared__20261107.npz',allow_pickle=False)
    ids=z['sample_index'];assert np.array_equal(ids,np.arange(len(s)))
    w=a['condition_weight'][ids]
    mae=mean_absolute_error(a['count72'][ids],z['count72_predicted'],sample_weight=w)
    mse=mean_squared_error(a['expression_z24'][ids],z['expression_predicted'],sample_weight=w)
    detail=pd.read_csv(OUT/'results/all_metrics.csv');row=detail[detail.protocol.eq('unseen_drug')&detail.method.eq('joint_shared')&detail.seed.astype(str).eq('20261107')&detail.scope.eq('all')].iloc[0]
    assert abs(mae-row.function_mae)<1e-8 and abs(mse-row.expression_mse)<1e-6
    summary=pd.read_csv(OUT/'results/seed_summary.csv');diff=pd.read_csv(OUT/'results/paired_drug_bootstrap.csv')
    for r in diff.itertuples():
        sub=summary[summary.protocol.eq(r.protocol)&summary.scope.eq('all')].set_index('method')
        expected_diff=sub.loc[r.method,r.metric+'_mean']-sub.loc[r.reference,r.metric+'_mean']
        assert abs(expected_diff-r.difference_mean)<1e-7
    ranks=pd.read_csv(OUT/'results/ranking_by_fold_cell_dose.csv')
    tied=ranks[ranks.method.isin(['train_mean','cell_mean','cell_dose_mean'])]
    assert np.allclose(tied.top5_observed_count,tied.all_drugs_mean_count,atol=1e-10)
    assert ranks.top5_regret.min()>-1e-10
    # Re-read original L5 values, including the additional 24h-only condition.
    prior=pd.read_csv(ROOT/'data/prepared_20261106/HMS_v3_实验准备/samples.csv')
    new=s[~s.library_condition_id.isin(prior.library_condition_id)]
    checks=sorted(set([0,len(s)-1]+new.index.tolist()))
    source=next((ROOT/'data/unpacked/LINCS_GSE70138').glob('*Level5*'))
    with h5py.File(source) as f:
        rows={v:i for i,v in enumerate(f['0/META/COL/id'].asstr()[:])}
        cols={v:i for i,v in enumerate(f['0/META/ROW/id'].asstr()[:])}
        gc=[cols[v] for v in a['gene_id']]
        for ix in checks:
            original=f['0/DATA/0/matrix'][rows[s.loc[ix,'sig_id']],:][gc]
            assert np.array_equal(original,a['expression_z24'][ix])
    source_audit=json.loads((OUT/'data/audit.json').read_text())
    for r in source_audit['source_hashes']:
        assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256']
    result={'neural_runs_checked':checkpoint_checks,'test_checkpoint_selection_violations':0,'condition_split_leaks':0,
            'sklearn_metric_crosscheck_passed':True,'bootstrap_differences_match_seed_means':True,
            'ranking_tie_expectation_check_passed':True,'dose_test_pairs_have_lower_training_doses':True,
            'source_L5_rows_reread':len(checks),'source_hashes_unchanged':True,
            'training_code_and_arrays_match_frozen_fingerprints':True,'new_training_in_verification':False}
    (OUT/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
