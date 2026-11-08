#!/usr/bin/env python3
"""Inspect existing data and replay ONE saved model; never train or alter data."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import torch
from train import OUT,ROOT,Model


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-index',type=int)
    args=parser.parse_args()
    samples=pd.read_csv(ROOT/'data/training/samples.csv')
    with np.load(ROOT/'data/training/arrays.npz',allow_pickle=False) as z:
        arrays={k:z[k] for k in z.files}
    assert np.array_equal(samples.library_condition_id.to_numpy(dtype=str),arrays['library_condition_id'])
    if args.sample_index is None:
        selected=samples[samples.cell_id.eq('BT20') & samples.pert_id.eq('BRD-A58767537') & np.isclose(samples.dose_uM,.04)]
        row=selected.iloc[0]
    else:
        row=samples[samples.sample_index.eq(args.sample_index)].iloc[0]
    i=int(row.sample_index)
    splits=pd.read_csv(ROOT/'data/training/splits.csv')
    roles=splits[splits.protocol.eq('unseen_drug') & splits.sample_index.eq(i)]
    fold=int(roles[roles.role.eq('test')].fold.iloc[0])
    method='joint_shared';seed=20261107
    run=OUT/'runs'/f'unseen_drug__fold{fold}__{method}__seed{seed}'
    with np.load(run/'preprocessing.npz',allow_pickle=False) as z:
        prep={k:z[k] for k in z.files}
    baseline=arrays['baseline'][int(row.cell_index)]
    drug=arrays['drug_morgan'][int(row.drug_index)]
    log_dose=float(arrays['log10_dose_uM'][i])
    cell_input=((baseline-prep['cell_mean'])/prep['cell_std']).astype(np.float32)
    dose_input=(log_dose-float(prep['dose_mean']))/float(prep['dose_std'])
    model=Model(method)
    model.load_state_dict(torch.load(run/'model_state.pt',map_location='cpu',weights_only=True))
    model.eval();torch.set_num_threads(2)
    with torch.no_grad():
        predicted_e,predicted_f,latent=model(torch.from_numpy(cell_input[None]),torch.tensor(drug[None],dtype=torch.float32),torch.tensor([[dose_input]],dtype=torch.float32))
    pe=predicted_e.numpy()[0]
    pf=float(predicted_f[0])*float(prep['function_std'])+float(prep['function_mean'])
    recorded=np.load(run/'test_predictions.npz',allow_pickle=False)
    j=int(np.flatnonzero(recorded['sample_index']==i)[0])
    gene_difference=float(np.max(abs(pe-recorded['expression_predicted'][j])))
    count_difference=abs(pf-float(recorded['count72_predicted'][j]))
    # Single-sample vs batched matrix multiplication can differ by float rounding.
    assert np.allclose(pe,recorded['expression_predicted'][j],atol=2e-5,rtol=2e-5)
    assert count_difference<2e-5
    export=ROOT/'examples'/f'sample_{i}'
    export.mkdir(parents=True,exist_ok=True)
    genes=pd.read_csv(ROOT/'data/raw/LINCS_GSE70138/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz',sep='\t').set_index('pr_gene_id')
    names=genes.loc[arrays['gene_id'].astype(int),'pr_gene_symbol'].to_numpy()
    pd.DataFrame({'gene_id':arrays['gene_id'],'gene_symbol':names,
                  'INPUT_baseline_log2TPMp1':baseline,'INPUT_baseline_after_train_scaling':cell_input,
                  'TARGET_observed_24h_response_z':arrays['expression_z24'][i],
                  'OUTPUT_predicted_24h_response_z':pe}).to_csv(export/'978基因_输入真实值与预测值.csv',index=False)
    pd.DataFrame({'bit_index':np.arange(len(drug)),'INPUT_drug_fingerprint':drug}).to_csv(export/'2048位药物指纹.csv',index=False)
    roles[['protocol','fold','role','heldout_group']].to_csv(export/'本样本在五折中的角色.csv',index=False)
    bio=pd.read_csv(ROOT/'data/reference/hms_biorepeats.csv')
    bio=bio[bio.library_condition_id.eq(row.library_condition_id)]
    assert abs(bio.relative_cell_count.mean()-float(arrays['count72'][i]))<1e-10
    bio.to_csv(export/'72h原始功能重复及孔引用.csv',index=False)
    (export/'本次模型结构.txt').write_text(str(model)+'\n')
    info={'sample_index':i,'cell':row.cell_id,'ModelID':row.ModelID,'drug':row.pert_iname,'pert_id':row.pert_id,
          'original_smiles':row.canonical_smiles,'library_id':row.library_id,'library_condition_id':row.library_condition_id,
          'cell_index_in_baseline_array':int(row.cell_index),'drug_index_in_fingerprint_array':int(row.drug_index),
          'dose_uM':float(row.dose_uM),'INPUT_log10_dose_uM':log_dose,'INPUT_dose_after_train_scaling':dose_input,
          'expression_signature_id':row.sig_id,'original_expression_instance_ids':row.distil_id,
          'expression_time_h':24,'function_time_h':72,'functional_label':'relative_cell_count_Nt_over_Ncontrol',
          'TARGET_observed_count72':float(arrays['count72'][i]),'OUTPUT_predicted_count72':pf,
          'replayed_run':str(run.relative_to(ROOT)),'method':method,'seed':seed,'fold':fold,'role_for_this_model':'test',
          'max_expression_difference_from_saved_test_prediction':gene_difference,
          'count_difference_from_saved_test_prediction':count_difference,
          'model_input_fields':['CCLE baseline scaled using this fold training statistics','Morgan fingerprint','log10 dose scaled using training statistics'],
          'target_fields_not_given_to_model':['observed 24h expression','observed 72h count'],
          'note':'This is one fold / one seed replay, not the three-seed ensemble. Only inspection exports were written; no optimizer or training was run.'}
    (export/'样本条件与前向核查.json').write_text(json.dumps(info,ensure_ascii=False,indent=2)+'\n')
    cells=samples[['cell_index','cell_id','ModelID']].drop_duplicates().sort_values('cell_index')
    compounds=samples[['drug_index','pert_id','pert_iname','canonical_smiles','chemical_split_group']].drop_duplicates().sort_values('drug_index')
    cells.to_csv(ROOT/'examples/细胞索引表.csv',index=False)
    compounds.to_csv(ROOT/'examples/药物索引表.csv',index=False)
    array_info=[{'array_key':k,'shape':str(v.shape),'dtype':str(v.dtype)} for k,v in arrays.items()]
    pd.DataFrame(array_info).to_csv(ROOT/'examples/训练数组形状.csv',index=False)
    print(json.dumps(info,ensure_ascii=False,indent=2))
    print('EXPORTED_TO',export)


if __name__=='__main__':main()
