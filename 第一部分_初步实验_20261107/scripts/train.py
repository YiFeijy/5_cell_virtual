#!/usr/bin/env python3
"""Frozen exploratory comparisons; independent run directories, resumable."""
import os
for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[variable] = '2'
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import hashlib
import json
import multiprocessing as mp
import random
import time
import numpy as np
import pandas as pd
import torch
from torch import nn

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '第一部分_初步实验_20261107'
METHODS = ['expr_only', 'func_only', 'joint_shared', 'joint_shared_E01', 'joint_shared_F01', 'joint_private', 'joint_bottleneck', 'function_bottleneck_control']
SEEDS = [20261107, 20261108, 20261109]


@lru_cache(maxsize=1)
def load():
    with np.load(OUT / 'data/arrays.npz', allow_pickle=False) as z:
        arrays = {k: z[k] for k in z.files}
    s = pd.read_csv(OUT / 'data/samples.csv')
    sp = pd.read_csv(OUT / 'data/splits.csv')
    assert np.array_equal(s.library_condition_id.to_numpy(dtype=str), arrays['library_condition_id'])
    return arrays, s, sp


@lru_cache(maxsize=6)
def inputs(protocol, fold):
    a, s, sp = load()
    d = sp[sp.protocol.eq(protocol) & sp.fold.eq(fold)].sort_values('sample_index')
    assert len(d) == len(s) and np.array_equal(d.library_condition_id.to_numpy(dtype=str), a['library_condition_id'])
    roles = {r: np.flatnonzero(d.role.eq(r)) for r in ['train', 'validation', 'test']}
    assert d.groupby('biological_condition_group').role.nunique().eq(1).all()
    if protocol == 'unseen_drug':
        assert d.groupby('heldout_group').role.nunique().eq(1).all()
    tr = roles['train']
    w = a['condition_weight']
    cell = a['baseline'][a['cell_index']].astype(np.float32)
    cm = np.average(cell[tr], axis=0, weights=w[tr])
    cs = np.sqrt(np.average((cell[tr]-cm)**2, axis=0, weights=w[tr])).clip(1e-4)
    dose = a['log10_dose_uM'][:, None]
    dm = float(np.average(dose[tr, 0], weights=w[tr]))
    ds = float(np.sqrt(np.average((dose[tr, 0]-dm)**2, weights=w[tr])))
    fm = float(np.average(a['count72'][tr], weights=w[tr]))
    fs = float(np.sqrt(np.average((a['count72'][tr]-fm)**2, weights=w[tr])))
    assert ds > 0 and fs > 0
    x = [((cell-cm)/cs).astype(np.float32), a['drug_morgan'][a['drug_index']].astype(np.float32), ((dose-dm)/ds).astype(np.float32)]
    assert all(np.isfinite(v).all() for v in x)
    preprocessing = {'cell_mean': cm, 'cell_std': cs, 'dose_mean': dm, 'dose_std': ds, 'function_mean': fm, 'function_std': fs}
    return x, a['expression_z24'].astype(np.float32), ((a['count72']-fm)/fs).astype(np.float32), w.astype(np.float32), roles, preprocessing


class Model(nn.Module):
    def __init__(self, method):
        super().__init__()
        self.method = method
        self.has_e = method not in ['func_only', 'function_bottleneck_control']
        self.has_f = method != 'expr_only'
        self.cell = nn.Sequential(nn.Linear(978, 128), nn.GELU())
        self.drug = nn.Sequential(nn.Linear(2048, 128), nn.GELU())
        self.dose = nn.Sequential(nn.Linear(1, 32), nn.GELU())
        self.fusion = nn.Sequential(nn.Linear(288, 128), nn.GELU(), nn.LayerNorm(128), nn.Dropout(.1))
        self.bottleneck = method in ['joint_bottleneck', 'function_bottleneck_control']
        if method == 'joint_private':
            self.shared = nn.Sequential(nn.Linear(128, 64), nn.GELU())
            self.private_e = nn.Sequential(nn.Linear(128, 64), nn.GELU())
            self.private_f = nn.Sequential(nn.Linear(128, 64), nn.GELU())
        self.expression = nn.Linear(128, 978) if self.has_e or self.bottleneck else None
        width = 978 if self.bottleneck else 128
        self.function = nn.Sequential(nn.Linear(width, 64), nn.GELU(), nn.Linear(64, 1)) if self.has_f else None

    def forward(self, cell, drug, dose):
        h = self.fusion(torch.cat([self.cell(cell), self.drug(drug), self.dose(dose)], -1))
        he = hf = h
        if self.method == 'joint_private':
            shared = self.shared(h)
            he = torch.cat([shared, self.private_e(h)], -1)
            hf = torch.cat([shared, self.private_f(h)], -1)
        e = self.expression(he) if self.expression is not None else None
        f = self.function(e if self.bottleneck else hf).squeeze(-1) if self.function is not None else None
        return e if self.has_e else None, f, e if self.bottleneck else h


def run(job):
    method, protocol, fold, seed = job
    path = OUT / 'runs' / f'{protocol}__fold{fold}__{method}__seed{seed}'
    if (path / 'complete.json').exists():
        return json.loads((path / 'complete.json').read_text())
    path.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); torch.set_num_threads(2)
    x, ye, yf, w, roles, prep = inputs(protocol, fold)
    x = [torch.from_numpy(v) for v in x]
    ye, yf, wt = torch.from_numpy(ye), torch.from_numpy(yf), torch.from_numpy(w)
    tr, va, te = (roles[k] for k in ['train', 'validation', 'test'])
    model = Model(method)
    opt = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=.0001)
    le_weight = .1 if method == 'joint_shared_E01' else 1.
    lf_weight = .1 if method == 'joint_shared_F01' else 1.
    def errors(ix, pe, pf):
        le = ((pe-ye[ix])**2).mean(1) if pe is not None else torch.zeros(len(ix))
        lf = (pf-yf[ix])**2 if pf is not None else torch.zeros(len(ix))
        return le, lf
    def evaluate(ix, output=False):
        model.eval()
        es, fs, zs, les, lfs = [], [], [], [], []
        with torch.no_grad():
            for start in range(0, len(ix), 512):
                ids = ix[start:start+512]
                pe, pf, latent = model(*(v[ids] for v in x))
                le, lf = errors(ids, pe, pf)
                les.append(le.numpy()); lfs.append(lf.numpy())
                if output:
                    if pe is not None: es.append(pe.numpy())
                    if pf is not None: fs.append(pf.numpy()*prep['function_std']+prep['function_mean'])
                    zs.append(latent.numpy())
        le, lf = np.concatenate(les), np.concatenate(lfs)
        # Same unweighted task-score definition for all joint weight variants.
        score = float(np.average(le+lf, weights=w[ix]))
        metrics = {'score': score, 'expression_mse': float(np.average(le, weights=w[ix])) if model.has_e else None,
                   'standardized_function_mse': float(np.average(lf, weights=w[ix])) if model.has_f else None}
        if output:
            return metrics, np.concatenate(es) if es else None, np.concatenate(fs) if fs else None, np.concatenate(zs)
        return metrics
    # A four-example forward/backward check precedes optimization.
    model.train(); ix = tr[:4]
    pe, pf, _ = model(*(v[ix] for v in x)); le, lf = errors(ix, pe, pf)
    initial_loss = ((le_weight*le+lf_weight*lf)*wt[ix]).sum()/wt[ix].sum()
    initial_loss.backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
    opt.zero_grad(set_to_none=True)
    best_score, best_epoch, best_state = float('inf'), 0, None
    history = []
    for epoch in range(1, 101):
        model.train(); order = np.random.permutation(tr); total, mass = 0., 0.
        for start in range(0, len(order), 256):
            ix = order[start:start+256]
            opt.zero_grad(set_to_none=True)
            pe, pf, _ = model(*(v[ix] for v in x)); le, lf = errors(ix, pe, pf)
            loss = ((le_weight*le+lf_weight*lf)*wt[ix]).sum()/wt[ix].sum()
            assert torch.isfinite(loss)
            loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 5.); opt.step()
            bw = float(wt[ix].sum()); total += float(loss.detach())*bw; mass += bw
        vm = evaluate(va)
        history.append({'epoch': epoch, 'train_weighted_loss': total/mass, **{f'validation_{k}': v for k, v in vm.items()}})
        if vm['score'] < best_score-1e-5:
            best_score, best_epoch = vm['score'], epoch
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        if epoch-best_epoch >= 15:
            break
    assert best_state is not None
    model.load_state_dict(best_state)
    for role, ix in [('validation', va), ('test', te)]:
        metrics, pe, pf, latent = evaluate(ix, output=True)
        payload = {'sample_index': ix, 'latent': latent}
        if pe is not None: payload['expression_predicted'] = pe
        if pf is not None: payload['count72_predicted'] = pf
        np.savez_compressed(path / f'{role}_predictions.npz', **payload)
        (path / f'{role}_normalized_loss.json').write_text(json.dumps(metrics, indent=2)+'\n')
    pd.DataFrame(history).to_csv(path / 'history.csv', index=False)
    torch.save(model.state_dict(), path / 'model_state.pt')
    np.savez_compressed(path / 'preprocessing.npz', **prep)
    result = {'method': method, 'protocol': protocol, 'fold': fold, 'seed': seed,
              'parameters': sum(p.numel() for p in model.parameters()), 'epochs_run': len(history), 'selected_epoch': best_epoch,
              'best_validation_score': best_score, 'expression_loss_weight': le_weight, 'function_loss_weight': lf_weight,
              'train_records': len(tr), 'validation_records': len(va), 'test_records': len(te),
              'seconds': time.monotonic()-started, 'device': 'cpu', 'finite_forward_backward_passed': True,
              'test_used_for_checkpoint_selection': False}
    (path / 'complete.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(); p.add_argument('--workers', type=int, default=6); p.add_argument('--verify-only', action='store_true'); a = p.parse_args()
    # Read all configuration and freeze the complete run list before any outer scoring.
    jobs = [(m, 'unseen_drug', f, seed) for m in METHODS for f in range(5) for seed in SEEDS]
    jobs += [(m, 'dose_extrapolation', 0, seed) for m in METHODS for seed in SEEDS]
    source = OUT / 'data/arrays.npz'
    frozen = {'date': '2026-11-07', 'main_task': '24h_L1000_z_and_72h_relative_cell_count', 'jobs': jobs,
              'max_epochs': 100, 'patience': 15, 'batch_size': 256, 'learning_rate': .001, 'weight_decay': .0001,
              'data_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'torch_version': torch.__version__}
    lock = OUT / 'frozen_training_protocol.json'
    if lock.exists():
        assert json.loads(lock.read_text()) == json.loads(json.dumps(frozen)), 'Frozen inputs/code/config changed; make a new experiment version.'
    else:
        lock.write_text(json.dumps(frozen, indent=2)+'\n')
    if a.verify_only:
        xx, ye, yf, w, roles, prep = inputs('unseen_drug', 0)
        for method in METHODS:
            model = Model(method); ix = roles['train'][:4]
            pe, pf, _ = model(*(torch.tensor(v[ix]) for v in xx))
            loss = torch.tensor(0.)
            if pe is not None: loss = loss+((pe-torch.tensor(ye[ix]))**2).mean()
            if pf is not None: loss = loss+((pf-torch.tensor(yf[ix]))**2).mean()
            loss.backward()
            assert all(v.grad is not None and torch.isfinite(v.grad).all() for v in model.parameters())
        print(json.dumps({'verified_architectures': METHODS, 'frozen_jobs': len(jobs), 'test_scored': False}), flush=True)
        return
    complete = []
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context('spawn')) as pool:
        futures = {pool.submit(run, job): job for job in jobs}
        for future in as_completed(futures):
            result = future.result(); complete.append(result)
            print(json.dumps({'completed': len(complete), 'total': len(jobs), **result}), flush=True)
            pd.DataFrame(complete).to_csv(OUT / 'training_progress.csv', index=False)
    pd.DataFrame(complete).sort_values(['protocol', 'method', 'fold', 'seed']).to_csv(OUT / 'training_runs.csv', index=False)


if __name__ == '__main__':
    main()
