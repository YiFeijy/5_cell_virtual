#!/usr/bin/env python3
"""Post-primary exploratory diagnostic with EXTRA measured 24h RNA input.

This is a different information budget, not an additional pre-treatment winner.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[v]='2'
import json
import numpy as np
import pandas as pd
from train import OUT,load,inputs
from analyze import metrics


def main():
    path=OUT/'information_diagnostic';path.mkdir(exist_ok=True)
    a,s,_=load();target=a['count72'];expr=a['expression_z24'].astype(np.float64)
    arrays={'measured_RNA_ridge':expr,'measured_RMS_ridge':np.sqrt((expr**2).mean(1))[:,None]}
    log=[];preds={key:np.full(len(s),np.nan) for key in arrays}
    for fold in range(5):
        _,_,_,w,roles,_=inputs('unseen_drug',fold);tr,va,te=(roles[r] for r in ['train','validation','test'])
        ym=float(np.average(target[tr],weights=w[tr]))
        for name,x in arrays.items():
            xm=np.average(x[tr],axis=0,weights=w[tr]);xc=x-xm
            gram=xc[tr].T@(w[tr,None]*xc[tr]);rhs=xc[tr].T@(w[tr]*(target[tr]-ym))
            eigen,u=np.linalg.eigh(gram);projected=u.T@rhs;best=None
            for alpha in [1.,10.,100.,1000.,10000.]:
                coef=u@(projected/(eigen+alpha))
                vp=xc[va]@coef+ym;score=float(np.average((vp-target[va])**2,weights=w[va]))
                log.append({'method':name,'fold':fold,'alpha':alpha,'validation_mse':score})
                if best is None or score<best[0]:best=score,alpha,coef
            preds[name][te]=xc[te]@best[2]+ym
            print(json.dumps({'method':name,'fold':fold,'selected_alpha':best[1]}),flush=True)
    rows=[]
    for name,pred in preds.items():
        assert np.isfinite(pred).all()
        np.savez_compressed(path/(name+'_OOF.npz'),sample_index=np.arange(len(s)),count72_predicted=pred)
        rows.append({'method':name,'input':'EXTRA measured 24h expression' if name=='measured_RNA_ridge' else 'EXTRA measured 24h response RMS',**metrics(np.arange(len(s)),None,pred)})
    pd.DataFrame(rows).to_csv(path/'metrics.csv',index=False);pd.DataFrame(log).to_csv(path/'validation_search.csv',index=False)
    (path/'scope.json').write_text(json.dumps({'status':'post-primary exploratory information diagnostic','primary_models_changed':False,'uses_extra_observation_at_prediction':True,'RNA_time_h':24,'function_time_h':72,'split':'same unseen_drug five-fold groups','comparison_interpretation':'Measures predictability given an actual post-treatment RNA assay. It is not fair to call it a superior pre-treatment model; association is not a causal state mechanism.'},indent=2)+'\n')
    print(pd.DataFrame(rows).to_string(index=False),flush=True)


if __name__=='__main__':main()
