#!/usr/bin/env python3
"""Validate the clean repository without training or reading large source matrices."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    count = 0
    manifest = ROOT / 'MANIFEST.csv'
    if manifest.exists():
        for r in csv.DictReader(manifest.open()):
            path = ROOT / r['path']
            assert path.is_file(), r['path']
            assert path.stat().st_size == int(r['bytes']), r['path']
            assert digest(path) == r['sha256'], r['path']
            count += 1
    frozen = json.loads((ROOT / 'results/frozen_training_protocol.json').read_text())
    assert digest(ROOT / 'src/train.py') == frozen['script_sha256']
    assert digest(ROOT / 'data/training/arrays.npz') == frozen['data_sha256']
    samples = pd.read_csv(ROOT / 'data/training/samples.csv')
    splits = pd.read_csv(ROOT / 'data/training/splits.csv')
    with np.load(ROOT / 'data/training/arrays.npz', allow_pickle=False) as z:
        assert len(samples) == 3364
        assert samples.biological_condition_group.nunique() == 3133
        assert np.array_equal(samples.sample_index, np.arange(len(samples)))
        assert np.array_equal(samples.library_condition_id.astype(str), z['library_condition_id'])
        assert z['baseline'].shape == (5,978) and z['drug_morgan'].shape == (106,2048)
        assert z['expression_z24'].shape == (3364,978)
        for k in ['baseline','expression_z24','count72','condition_weight','log10_dose_uM']:
            assert np.isfinite(z[k]).all(), k
        assert np.allclose(samples.groupby('biological_condition_group').condition_weight.sum(), 1)
    for (protocol, fold), group in splits.groupby(['protocol','fold']):
        assert set(group.sample_index) == set(samples.sample_index)
        assert group.groupby('biological_condition_group').role.nunique().eq(1).all()
        if protocol in ['unseen_drug','unseen_scaffold']:
            assert group.groupby('heldout_group').role.nunique().eq(1).all()
    runs = pd.read_csv(ROOT / 'results/training_runs.csv')
    history = pd.read_csv(ROOT / 'results/training_history.csv.gz')
    assert len(runs) == history.run_id.nunique() == 144
    assert len(history) == 4724
    assert not runs.test_used_for_checkpoint_selection.any()
    predictions = pd.read_csv(ROOT / 'results/metrics/functional_predictions.csv.gz', dtype={'seed':str})
    metrics = pd.read_csv(ROOT / 'results/metrics/all_metrics.csv', dtype={'seed':str})
    groups = 0
    for (protocol,method,seed), group in predictions.groupby(['protocol','method','seed']):
        row = metrics[metrics.protocol.eq(protocol)&metrics.method.eq(method)&metrics.seed.eq(seed)&metrics.scope.eq('all')]
        assert len(row) == 1
        error = np.average(abs(group.observed_count72-group.predicted_count72), weights=group.condition_weight)
        assert abs(error-row.function_mae.iloc[0]) < 1e-10
        groups += 1
    print(json.dumps({'file_hashes_checked':count,'observations':len(samples),'conditions':3133,
                      'split_checks':'passed','historical_runs':len(runs),'epochs':len(history),
                      'functional_prediction_groups_verified':groups,'training_executed':False}, indent=2))


if __name__ == '__main__':
    main()
