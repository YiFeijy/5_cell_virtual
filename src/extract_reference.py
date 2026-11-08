#!/usr/bin/env python3
"""Extract the fixed 24h candidate signatures from the official Level5 matrix."""
from pathlib import Path
import json
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    reference = ROOT/'data/reference'
    signatures = pd.read_csv(reference/'lincs_signatures.csv')
    with np.load(reference/'input_features.npz', allow_pickle=False) as z:
        genes = z['gene_id']
    annotation = pd.read_csv(next((ROOT/'data/raw/LINCS_GSE70138').glob('*gene_info*')), sep='\t')
    assert set(genes.astype(int)) == set(annotation.loc[annotation.pr_is_lm.eq(1),'pr_gene_id'])
    source = next((ROOT/'data/unpacked/LINCS_GSE70138').glob('*Level5*.gctx'))
    result = np.empty((len(signatures),len(genes)),dtype=np.float32)
    with h5py.File(source,'r') as f:
        row_ids = f['0/META/COL/id'].asstr()[:]
        column_ids = f['0/META/ROW/id'].asstr()[:]
        row_map = {v:i for i,v in enumerate(row_ids)}
        col_map = {v:i for i,v in enumerate(column_ids)}
        cols = np.array([col_map[g] for g in genes])
        ordering = sorted(range(len(signatures)),key=lambda i:row_map[signatures.sig_id.iloc[i]])
        for start in range(0,len(ordering),256):
            selected = ordering[start:start+256]
            block = f['0/DATA/0/matrix'][[row_map[signatures.sig_id.iloc[i]] for i in selected],:]
            result[selected] = block[:,cols]
    target = reference/'lincs_expression.npz'
    matched_existing = target.exists()
    if matched_existing:
        with np.load(target,allow_pickle=False) as previous:
            assert np.array_equal(previous['sig_id'],signatures.sig_id.to_numpy(dtype=str))
            assert np.array_equal(previous['expression_z'],result)
    assert np.isfinite(result).all()
    np.savez_compressed(target,expression_z=result,sig_id=signatures.sig_id.to_numpy(dtype=str),gene_id=genes)
    print(json.dumps({'source':source.name,'signatures':len(signatures),'measured_genes':len(genes),'existing_values_match':matched_existing}))


if __name__ == '__main__':
    main()
