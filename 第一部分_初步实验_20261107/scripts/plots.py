#!/usr/bin/env python3
"""Standalone scientific figures for exploratory first-part results."""
from pathlib import Path
import os
os.environ['MPLCONFIGDIR'] = '/tmp/mvp_matplotlib_20261107'
os.environ['XDG_CACHE_HOME'] = '/tmp/mvp_cache_20261107'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from train import OUT, load

LABELS = {'zero_response':'Zero response', 'train_mean':'Global mean', 'cell_mean':'Cell mean', 'cell_dose_mean':'Cell + dose mean',
          'ridge_independent':'Ridge', 'expr_only':'ST expression', 'func_only':'ST function', 'joint_shared':'Shared MTL',
          'joint_shared_E01':'Shared MTL (E x0.1)', 'joint_shared_F01':'Shared MTL (F x0.1)', 'joint_private':'Shared / private',
          'joint_bottleneck':'Expression bottleneck', 'function_bottleneck_control':'Bottleneck capacity control'}
COLORS = {'func_only':'#2563a6', 'expr_only':'#2563a6', 'joint_shared':'#df7d1e', 'joint_bottleneck':'#26856d', 'joint_private':'#7b5ca9'}


def save(fig, name):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/(name+'.png'),dpi=180,bbox_inches='tight')
    fig.savefig(folder/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)


def prediction(protocol, method):
    suffix='deterministic' if method in ['zero_response','train_mean','cell_mean','cell_dose_mean','ridge_independent'] else 'ensemble'
    return np.load(OUT/'predictions'/f'{protocol}__{method}__{suffix}.npz',allow_pickle=False)


def comparison(summary, protocol, filename):
    d=summary[summary.protocol.eq(protocol)&summary.scope.eq('all')].set_index('method')
    fig,axes=plt.subplots(1,2,figsize=(12,6.6))
    for ax,metric,title in zip(axes,['function_mae','expression_mse'],['72h relative cell count: MAE','24h response: MSE']):
        names=[k for k in LABELS if k in d.index and np.isfinite(d.loc[k,metric+'_mean'])]
        values=d.loc[names,metric+'_mean'];err=d.loc[names,metric+'_std'].fillna(0)
        ax.barh(range(len(names)),values,xerr=err,color=[COLORS.get(k,'#91a0ad') for k in names],alpha=.9,capsize=2)
        ax.set_yticks(range(len(names)),[LABELS[k] for k in names],fontsize=8.5);ax.invert_yaxis();ax.set_xlabel('Error (lower is better)');ax.set_title(title,fontsize=11)
        ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    fig.suptitle('Unseen-drug 5-fold OOF' if protocol=='unseen_drug' else 'Dose extrapolation: train <=1.11 uM; test 10 uM',fontsize=13)
    fig.text(.5,.01,'Neural bars: mean +/- SD across 3 seeds; deterministic baselines: one result. Current five-cell cohort.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.03,1,.96]);save(fig,filename)


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
    a,s,_=load();r=OUT/'results';summary=pd.read_csv(r/'seed_summary.csv')
    comparison(summary,'unseen_drug','01_unseen_drug_comparison')
    comparison(summary,'dose_extrapolation','02_dose_extrapolation_comparison')
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,method in zip(axes,['func_only','joint_shared','joint_bottleneck']):
        p=prediction('unseen_drug',method);ix=p['sample_index'];truth=a['count72'][ix];pred=p['count72_predicted']
        ax.scatter(truth,pred,s=4,alpha=.22,rasterized=True,color=COLORS[method]);low=min(-.05,float(pred.min()));high=max(1.25,float(pred.max()))
        ax.plot([low,high],[low,high],'k--',lw=1);ax.set(xlim=(low,high),ylim=(low,high),xlabel='Observed relative count',ylabel='Predicted relative count',title=LABELS[method]);ax.set_aspect('equal')
    fig.suptitle('Unseen-drug OOF predictions; 3-seed ensembles',fontsize=12);fig.tight_layout();save(fig,'03_function_predictions')
    fig,axes=plt.subplots(1,2,figsize=(11,3.7))
    for ax,metric,methods in [(axes[0],'validation_expression_mse',['expr_only','joint_shared','joint_private','joint_bottleneck']),
                              (axes[1],'validation_standardized_function_mse',['func_only','joint_shared','joint_private','joint_bottleneck'])]:
        for method in methods:
            p=OUT/'runs'/f'unseen_drug__fold0__{method}__seed20261107'/'history.csv';h=pd.read_csv(p)
            ax.plot(h.epoch,h[metric],label=LABELS[method],color=COLORS[method],lw=1.5)
        ax.set_xlabel('Epoch');ax.set_ylabel('Validation MSE');ax.set_title('24h response' if ax==axes[0] else '72h count, standardized training scale');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Training diagnostic: fold 0, seed 20261107; no test-based epoch selection',fontsize=11);fig.tight_layout();save(fig,'04_validation_curves')
    # Examples fixed by earlier audited records, not selected for low prediction error.
    fig,axes=plt.subplots(1,2,figsize=(10,3.9))
    for ax,name in zip(axes,['afatinib','GSK-1059615']):
        d=s[s.cell_id.eq('BT20') & s.pert_iname.str.lower().eq(name.lower())].copy()
        if d.empty:continue
        observed=d.groupby('dose_uM').relative_cell_count.mean()
        ax.plot(observed.index,observed,'ko-',lw=1.2,label='Observed library mean')
        for method in ['func_only','joint_shared','joint_bottleneck']:
            p=prediction('unseen_drug',method);lookup=pd.Series(p['count72_predicted'],index=p['sample_index'])
            dd=d.copy();dd['pred']=lookup.reindex(dd.sample_index).to_numpy();curve=dd.groupby('dose_uM').pred.mean()
            ax.plot(curve.index,curve,'o--',ms=3,lw=1.2,label=LABELS[method],color=COLORS[method])
        ax.set_xscale('log');ax.set_xlabel('Dose (uM)');ax.set_ylabel('72h relative cell count');ax.set_title(f'BT20 + {name}');ax.legend(fontsize=7.5);ax.grid(alpha=.2)
    fig.suptitle('Preselected worked examples; unseen-drug OOF ensembles',fontsize=12);fig.tight_layout();save(fig,'05_dose_curves_examples')
    ranks=pd.read_csv(r/'ranking_seed_summary.csv');ranks=ranks[ranks.protocol.eq('unseen_drug')&ranks.seed.astype(str).ne('ensemble')]
    rr=ranks.groupby('method').top5_regret.agg(['mean','std'])
    names=[k for k in LABELS if k in rr.index]
    fig,ax=plt.subplots(figsize=(8,5));ax.barh(range(len(names)),rr.loc[names,'mean'],xerr=rr.loc[names,'std'].fillna(0),color=[COLORS.get(k,'#91a0ad') for k in names],capsize=2)
    ax.set_yticks(range(len(names)),[LABELS[k] for k in names],fontsize=9);ax.invert_yaxis();ax.set_xlabel('Observed count regret of selected top 5 (lower is better)');ax.set_title('Screening within each held-out drug fold, cell and dose');ax.grid(axis='x',alpha=.2);fig.tight_layout();save(fig,'06_top5_screening_regret')
    differences=pd.read_csv(r/'paired_drug_bootstrap.csv');d=differences[differences.protocol.eq('unseen_drug')&differences.metric.eq('function_mae')&differences.reference.eq('func_only')]
    fig,ax=plt.subplots(figsize=(8,3.5));y=np.arange(len(d));ax.errorbar(d.difference_mean,y,xerr=np.array([d.difference_mean-d.ci95_low,d.ci95_high-d.difference_mean]),fmt='o',capsize=3,color='#385a73');ax.axvline(0,color='k',ls='--',lw=1)
    ax.set_yticks(y,[LABELS[k] for k in d.method]);ax.invert_yaxis();ax.set_xlabel('MAE difference vs ST function; negative favors joint model');ax.set_title('Paired drug bootstrap: 95% descriptive intervals');ax.grid(axis='x',alpha=.2);fig.tight_layout();save(fig,'07_paired_function_differences')
    diagnostic=OUT/'information_diagnostic/metrics.csv'
    if diagnostic.exists():
        diag=pd.read_csv(diagnostic).set_index('method')
        values=[float(summary[summary.protocol.eq('unseen_drug')&summary.scope.eq('all')&summary.method.eq('func_only')].function_mae_mean.iloc[0]),
                float(diag.loc['measured_RMS_ridge','function_mae']),float(diag.loc['measured_RNA_ridge','function_mae'])]
        fig,ax=plt.subplots(figsize=(7,3.5));ax.barh(range(3),values,color=['#2563a6','#91a0ad','#26856d'])
        ax.set_yticks(range(3),['Pre-treatment inputs: ST function','EXTRA measured 24h response magnitude','EXTRA measured 24h full RNA signature']);ax.invert_yaxis();ax.set_xlabel('72h count MAE on held-out drugs');ax.set_title('Information diagnostic: different measurement budgets')
        for i,v in enumerate(values):ax.text(v+.002,i,f'{v:.4f}',va='center',fontsize=9)
        fig.tight_layout();save(fig,'08_measured_RNA_information_diagnostic')
    print('Saved main figures and available diagnostic, each in PNG and PDF.')


if __name__=='__main__':
    main()
