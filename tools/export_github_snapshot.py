#!/usr/bin/env python3
"""Copy a curated research snapshot to an existing, separate Git checkout.

No git add/commit/push and no deletions are performed by this script.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = '第一部分_初步实验_20261107'
SELECTED = 'unseen_drug__fold1__joint_shared__seed20261107'
PREDICTION = 'unseen_drug__joint_shared__20261107.npz'


def include(path):
    parts = path.parts
    if any(p.startswith('.') and p not in {'.gitignore'} for p in parts):
        return False
    if '__pycache__' in parts or path.suffix in {'.pyc', '.part', '.tmp', '.lock'}:
        return False
    if path.name.startswith('github_发布回执_'):
        return False
    if parts[:2] in [('data', 'unpacked'), ('data', 'download_logs')]:
        return False
    if parts[:2] == ('data', 'raw'):
        dataset = parts[2]
        if dataset == 'HMS_LINCS':
            return path.suffix == '.xls'
        if dataset == 'LINCS_GSE70138':
            return path.name.endswith('.txt.gz')
        if dataset == 'CCLE_DepMap_24Q4':
            return path.name in {'Model.csv', 'OmicsProfiles.csv', 'OmicsDefaultModelProfiles.csv'}
        return False
    if parts[:2] == ('数据可行性调研_20261105', 'metadata'):
        return path.suffix not in {'.zip', '.rds'} and (ROOT / path).stat().st_size < 12 * 1024**2
    if parts[:2] == (EXPERIMENT, 'runs'):
        return parts[2] == SELECTED
    if parts[:2] == (EXPERIMENT, 'baselines'):
        return False
    if parts[:2] == (EXPERIMENT, 'predictions'):
        return path.name == PREDICTION
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    target = args.destination.resolve()
    assert target != ROOT and not target.is_relative_to(ROOT), 'Use a separate checkout outside the source project.'
    assert (target / '.git').exists(), 'Destination must already be a Git checkout.'
    copied = []
    for source in sorted(ROOT.rglob('*')):
        if not source.is_file() or source.is_symlink():
            continue
        relative = source.relative_to(ROOT)
        if not include(relative):
            continue
        assert source.stat().st_size < 50 * 1024**2, f'Review this large file explicitly: {relative}'
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        copied.append(relative)
    # Historical complete.json files are kept away from resumable runs/.
    records = target / EXPERIMENT / 'run_records'
    for run in sorted((ROOT / EXPERIMENT / 'runs').iterdir()):
        if not run.is_dir():
            continue
        for source in sorted(run.iterdir()):
            if source.name in {'complete.json', 'history.csv', 'preprocessing.npz',
                               'test_normalized_loss.json', 'validation_normalized_loss.json'}:
                destination = records / run.name / source.name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
    samples = pd.read_csv(ROOT / EXPERIMENT / 'data/samples.csv')
    with np.load(ROOT / EXPERIMENT / 'data/arrays.npz', allow_pickle=False) as arrays:
        true_count = arrays['count72']
    predictions = []
    for source in sorted((ROOT / EXPERIMENT / 'predictions').glob('*.npz')):
        protocol, method, seed = source.stem.split('__')
        with np.load(source, allow_pickle=False) as arrays:
            if 'count72_predicted' not in arrays.files:
                continue
            indices = arrays['sample_index']
            frame = samples.iloc[indices][['sample_index', 'library_condition_id', 'biological_condition_group',
                                           'cell_id', 'pert_id', 'pert_iname', 'dose_uM', 'condition_weight']].copy()
            frame['protocol'], frame['method'], frame['seed'] = protocol, method, seed
            frame['observed_count72'] = true_count[indices]
            frame['predicted_count72'] = arrays['count72_predicted']
            predictions.append(frame)
    joined = pd.concat(predictions, ignore_index=True)
    joined.to_csv(target / EXPERIMENT / 'results/functional_predictions.csv.gz', index=False,
                  compression={'method': 'gzip', 'mtime': 0})
    manifest = target / 'GITHUB_CONTENTS.csv'
    with manifest.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['path', 'bytes', 'sha256'])
        for path in sorted(target.rglob('*')):
            if not path.is_file() or '.git' in path.parts or '__pycache__' in path.parts or path == manifest:
                continue
            digest = hashlib.sha256()
            with path.open('rb') as reader:
                for block in iter(lambda: reader.read(1024 * 1024), b''):
                    digest.update(block)
            writer.writerow([str(path.relative_to(target)), path.stat().st_size, digest.hexdigest()])
    print(json.dumps({'destination': str(target), 'copied_files': len(copied),
                      'historical_run_records': len(list(records.iterdir())),
                      'functional_prediction_rows': len(joined), 'manifest': str(manifest)},
                     ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
