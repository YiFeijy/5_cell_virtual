#!/usr/bin/env python3
"""One pre-QC train/validation pilot, without scoring the outer test set.

This checks whether the framework's first-part joint dataflow can run. It is
not a completed comparison, a final training release or a survivor-state model.
"""
from pathlib import Path
import argparse
import hashlib
import json
import random
import time
import numpy as np
import pandas as pd
import torch
from torch import nn

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/prepared_20261106/HMS_v3_实验准备'
OUT = ROOT / '第一部分_依据整体框架的设计与起步_20261106/联合模型调通试验'


class JointStateModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.cell = nn.Sequential(nn.Linear(978, 128), nn.GELU())
        self.drug = nn.Sequential(nn.Linear(2048, 128), nn.GELU())
        self.dose = nn.Sequential(nn.Linear(1, 32), nn.GELU())
        self.context = nn.Sequential(nn.Linear(288, 128), nn.GELU())
        self.time = nn.Sequential(nn.Linear(1, 16), nn.GELU())
        self.state = nn.Sequential(nn.Linear(144, 128), nn.GELU(), nn.LayerNorm(128))
        self.expression = nn.Linear(128, 978)
        self.function72 = nn.Sequential(nn.Linear(256, 64), nn.GELU(), nn.Linear(64, 1))
        self.register_buffer('expression_queries', torch.log1p(torch.tensor([3., 24.])) / np.log(25.))

    def forward(self, cell, drug, dose):
        context = self.context(torch.cat([self.cell(cell), self.drug(drug), self.dose(dose)], dim=-1))
        tq = self.time(self.expression_queries[:, None])
        state = self.state(torch.cat([context[:, None].expand(-1, 2, -1), tq[None].expand(len(cell), -1, -1)], dim=-1))
        return self.expression(state), self.function72(state.flatten(1)).squeeze(-1), state


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epochs', type=int, default=12)
    p.add_argument('--seed', type=int, default=20261106)
    a = p.parse_args()
    assert a.epochs > 0
    torch.set_num_threads(2)
    random.seed(a.seed)
    np.random.seed(a.seed)
    torch.manual_seed(a.seed)
    OUT.mkdir(exist_ok=True, parents=True)
    arrays = np.load(DATA / 'arrays.npz', allow_pickle=False)
    samples = pd.read_csv(DATA / 'samples.csv')
    splits = pd.read_csv(DATA / 'split_assignments_DRAFT.csv')
    splits = splits[splits.protocol.eq('unseen_drug') & splits.fold.eq(0)].sort_values('sample_index')
    assert np.array_equal(samples.library_condition_id.to_numpy(dtype=str), arrays['library_condition_id'])
    assert np.array_equal(splits.library_condition_id.to_numpy(dtype=str), arrays['library_condition_id'])
    assert splits.groupby('biological_condition_group').role.nunique().eq(1).all()
    train = np.flatnonzero(splits.role.eq('train'))
    val = np.flatnonzero(splits.role.eq('validation'))
    test = np.flatnonzero(splits.role.eq('test'))
    assert set(samples.iloc[train].chemical_split_group).isdisjoint(samples.iloc[val].chemical_split_group)
    assert set(samples.iloc[train].chemical_split_group).isdisjoint(samples.iloc[test].chemical_split_group)
    weights = arrays['condition_weight']
    cell = arrays['baseline_log2TPMp1'][arrays['cell_index']].astype(np.float32)
    drug = arrays['drug_morgan'][arrays['drug_index']].astype(np.float32)
    dose = arrays['log10_dose_uM'][:, None].astype(np.float32)
    mean = np.average(cell[train], axis=0, weights=weights[train])
    std = np.sqrt(np.average((cell[train] - mean) ** 2, axis=0, weights=weights[train])).clip(1e-4)
    dose_mean = float(np.average(dose[train, 0], weights=weights[train]))
    dose_std = float(np.sqrt(np.average((dose[train, 0] - dose_mean) ** 2, weights=weights[train])))
    func_mean = float(np.average(arrays['relative_cell_count_72h'][train], weights=weights[train]))
    func_std = float(np.sqrt(np.average((arrays['relative_cell_count_72h'][train] - func_mean) ** 2, weights=weights[train])))
    assert dose_std > 0 and func_std > 0
    features = [torch.from_numpy(((cell - mean) / std).astype(np.float32)), torch.from_numpy(drug), torch.from_numpy(((dose - dose_mean) / dose_std).astype(np.float32))]
    expression = torch.from_numpy(arrays['expression_z'].astype(np.float32))
    function = torch.from_numpy(((arrays['relative_cell_count_72h'] - func_mean) / func_std).astype(np.float32))
    wt = torch.from_numpy(weights.astype(np.float32))
    model = JointStateModel()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)

    def loss_for(ix, pred_e, pred_f):
        le = ((pred_e - expression[ix]) ** 2).mean((1, 2))
        lf = (pred_f - function[ix]) ** 2
        return ((le + lf) * wt[ix]).sum() / wt[ix].sum()

    def evaluate(ix):
        model.eval()
        with torch.no_grad():
            pe, pf, state = model(*(x[ix] for x in features))
            score = float(loss_for(ix, pe, pf))
            emse = float(np.average(((pe.numpy() - expression[ix].numpy()) ** 2).mean((1, 2)), weights=weights[ix]))
            raw_f = pf.numpy() * func_std + func_mean
            fmae = float(np.average(abs(raw_f - arrays['relative_cell_count_72h'][ix]), weights=weights[ix]))
        return {'loss': score, 'expression_mse': emse, 'count72_mae': fmae}, pe.numpy(), raw_f, state.numpy()

    # Explicit finite-gradient check for the two supervised outputs.
    model.train()
    ix = train[:4]
    pe, pf, state = model(*(x[ix] for x in features))
    assert pe.shape == (4, 2, 978) and pf.shape == (4,) and state.shape == (4, 2, 128)
    loss_for(ix, pe, pf).backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
    optimizer.zero_grad(set_to_none=True)
    initial, _, _, _ = evaluate(val)
    history = [{'epoch': 0, 'train_loss': None, **{f'validation_{k}': v for k, v in initial.items()}}]
    best = (initial['loss'], 0, {k: v.detach().clone() for k, v in model.state_dict().items()})
    started = time.monotonic()
    for epoch in range(1, a.epochs + 1):
        model.train()
        total, total_w = 0., 0.
        order = np.random.permutation(train)
        for start in range(0, len(order), 128):
            ix = order[start:start + 128]
            optimizer.zero_grad(set_to_none=True)
            pe, pf, _ = model(*(x[ix] for x in features))
            loss = loss_for(ix, pe, pf)
            assert torch.isfinite(loss)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.)
            optimizer.step()
            batch_w = float(wt[ix].sum())
            total += float(loss.detach()) * batch_w
            total_w += batch_w
        score, _, _, _ = evaluate(val)
        history.append({'epoch': epoch, 'train_loss': total / total_w, **{f'validation_{k}': v for k, v in score.items()}})
        if score['loss'] < best[0]:
            best = (score['loss'], epoch, {k: v.detach().clone() for k, v in model.state_dict().items()})
        print(json.dumps({'epoch': epoch, 'train_loss': total / total_w, 'validation_loss': score['loss']}), flush=True)
    model.load_state_dict(best[2])
    final, pe, pf, state = evaluate(val)
    pd.DataFrame(history).to_csv(OUT / 'training_history.csv', index=False)
    torch.save(model.state_dict(), OUT / 'model_state.pt')
    np.savez_compressed(OUT / 'preprocessing.npz', baseline_mean=mean, baseline_std=std,
                        dose_mean=dose_mean, dose_std=dose_std, function_mean=func_mean, function_std=func_std,
                        gene_id=arrays['gene_id'])
    np.savez_compressed(OUT / 'validation_predictions.npz', expression_predicted=pe, function_predicted=pf,
                        latent_response_state=state, library_condition_id=arrays['library_condition_id'][val])
    vp = samples.iloc[val][['library_condition_id', 'biological_condition_group', 'cell_id', 'pert_id', 'dose_uM', 'condition_weight', 'relative_cell_count']].copy()
    vp['predicted_count72'] = pf
    vp.to_csv(OUT / 'validation_predictions.csv', index=False)
    result = {'stage': 'PRE_QC_DATAFLOW_PILOT_train_validation_only', 'model_fits': 1, 'seed': a.seed,
              'epochs': a.epochs, 'selected_epoch_by_validation': best[1], 'device': 'cpu',
              'train_library_records': len(train), 'validation_library_records': len(val),
              'outer_test_library_records_untouched': len(test), 'outer_test_scored': False,
              'parameter_count': sum(p.numel() for p in model.parameters()),
              'expression_output_shape_per_condition': [2, 978], 'state_shape_per_condition': [2, 128],
              'functional_output': '72h relative cell count, one value per library condition',
              'initial_validation': initial, 'selected_validation': final,
              'finite_forward_and_backward_passed': True, 'feature_scaling_fit_on_train_only': True,
              'observed_post_treatment_expression_used_as_input': False,
              'observed_function_used_as_input': False, 'biological_QC_finalized': False,
              'independent_survivor_state_reference_used': False, 'baseline_comparison_completed': False,
              'seconds': time.monotonic() - started, 'torch_version': torch.__version__,
              'data_audit_source': str((DATA / 'audit.json').relative_to(ROOT)),
              'sample_manifest_sha256': hashlib.sha256((DATA / 'samples.csv').read_bytes()).hexdigest(),
              'interpretation': 'A runnable shared response-state model only. Validation was used for checkpoint selection. This run is not an outer-test result, a biological state validation or evidence of superiority/novelty.'}
    (OUT / 'run_summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


if __name__ == '__main__':
    main()
