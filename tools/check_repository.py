#!/usr/bin/env python3
"""Check the published training data, split integrity and representative run."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / '第一部分_初步实验_20261107'


def sha256(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    manifest = ROOT / 'GITHUB_CONTENTS.csv'
    verified = 0
    if manifest.exists():
        with manifest.open() as stream:
            for row in csv.DictReader(stream):
                path = ROOT / row['path']
                assert path.is_file(), f'Missing published file: {row["path"]}'
                assert path.stat().st_size == int(row['bytes']), row['path']
                assert sha256(path) == row['sha256'], row['path']
                verified += 1
    frozen = json.loads((EXP / 'frozen_training_protocol.json').read_text())
    assert sha256(EXP / 'data/arrays.npz') == frozen['data_sha256']
    assert sha256(EXP / 'scripts/train.py') == frozen['script_sha256']
    samples = pd.read_csv(EXP / 'data/samples.csv')
    splits = pd.read_csv(EXP / 'data/splits.csv')
    with np.load(EXP / 'data/arrays.npz', allow_pickle=False) as arrays:
        assert len(samples) == 3364
        assert np.array_equal(samples.sample_index, np.arange(len(samples)))
        assert np.array_equal(samples.library_condition_id.astype(str), arrays['library_condition_id'])
        assert arrays['expression_z24'].shape == (3364, 978)
        assert arrays['baseline'].shape == (5, 978)
        assert arrays['drug_morgan'].shape == (106, 2048)
        assert np.array_equal(samples.cell_index, arrays['cell_index'])
        assert np.array_equal(samples.drug_index, arrays['drug_index'])
        for key in ['expression_z24', 'baseline', 'count72', 'log10_dose_uM', 'condition_weight']:
            assert np.isfinite(arrays[key]).all(), key
        sums = pd.DataFrame({'condition': samples.biological_condition_group,
                             'weight': arrays['condition_weight']}).groupby('condition').weight.sum()
        assert len(sums) == 3133 and np.allclose(sums, 1)
        for (protocol, fold), group in splits.groupby(['protocol', 'fold']):
            assert set(group.sample_index) == set(samples.sample_index)
            assert group.groupby('biological_condition_group').role.nunique().eq(1).all()
            if protocol in ['unseen_drug', 'unseen_scaffold']:
                assert group.groupby('heldout_group').role.nunique().eq(1).all()
        representative = EXP / 'runs/unseen_drug__fold1__joint_shared__seed20261107'
        for name in ['model_state.pt', 'preprocessing.npz', 'test_predictions.npz', 'history.csv']:
            assert (representative / name).is_file(), name
        with np.load(representative / 'test_predictions.npz', allow_pickle=False) as prediction:
            expected = splits.query('protocol == "unseen_drug" and fold == 1 and role == "test"')
            assert set(prediction['sample_index']) == set(expected.sample_index)
        runs = pd.read_csv(EXP / 'training_runs.csv')
        assert len(runs) == 144
        assert not runs.test_used_for_checkpoint_selection.any()
        if (EXP / 'run_records').exists():
            assert len(list((EXP / 'run_records').glob('*/history.csv'))) == 144
        result = {'published_file_hashes_checked': verified, 'observations': len(samples),
                  'distinct_conditions': len(sums), 'cells': 5, 'drugs': 106,
                  'split_integrity_passed': True, 'frozen_data_and_train_code_match': True,
                  'representative_checkpoint_present': True, 'historical_training_runs': len(runs),
                  'training_executed': False}
        print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
