"""Frozen external evaluation. No fitting or model selection. Run from this directory."""
from pathlib import Path
import json, hashlib, zipfile, io, shutil
import numpy as np
import pandas as pd
import joblib
from scipy.special import expit, softmax
from sklearn.metrics import roc_auc_score, average_precision_score, log_loss
ROOT=Path(__file__).resolve().parent
REPS=1000

def calibrated(m,p):
 c=m['calibration']; v=c['parameters']
 if c['method']=='Platt_sigmoid_validation':
  q=np.clip(p[:,list(m['classes']).index(1)],1e-7,1-1e-7)
  q=expit(np.log(q/(1-q))*float(v['coef'][0][0])+float(v['intercept'][0]))
  return np.column_stack([1-q,q])
 if c['method']=='temperature_validation':return softmax(np.log(np.clip(p,1e-7,1))/float(v['temperature']),axis=1)
 raise ValueError(c)

def metrics(y,p):
 k=p.shape[1];oh=np.eye(k)[y]
 if k==2:
  auc=roc_auc_score(y,p[:,1]);ap=average_precision_score(y,p[:,1]);bs=np.mean((y-p[:,1])**2);conf=p[:,1];obs=y
 else:
  auc=roc_auc_score(y,p,multi_class='ovr',average='macro',labels=np.arange(k));ap=average_precision_score(oh,p,average='macro');bs=np.mean(np.sum((oh-p)**2,axis=1));conf=p.max(1);obs=(p.argmax(1)==y).astype(float)
 bins=np.minimum((conf*10).astype(int),9);ece=sum(np.mean(bins==b)*abs(np.mean(conf[bins==b])-np.mean(obs[bins==b])) for b in range(10) if np.any(bins==b))
 return dict(AUROC=float(auc),AP=float(ap),Brier=float(bs),LogLoss=float(log_loss(y,p,labels=np.arange(k))),ECE10=float(ece))

def main():
 x=np.load(ROOT/'inputs/patient_features.npy');idx=pd.read_csv(ROOT/'inputs/patient_index.csv');lab=pd.read_csv(ROOT/'inputs/endpoint_reference_freeze.csv')
 assert x.shape==(len(idx),1024) and np.isfinite(x).all() and idx.patient_id.is_unique
 assert not lab.duplicated(['patient_id','endpoint']).any()
 labels=lab[lab.patient_id.isin(idx.patient_id)].copy();assert labels.patient_id.nunique()==len(idx)
 labels.to_csv(ROOT/'external_label_alignment.csv',index=False)
 rows=[];preds=[];elig=[];bins_out=[]
 for ep in sorted(labels.endpoint.unique()):
  ll=labels[labels.endpoint.eq(ep)].set_index('patient_id').reindex(idx.patient_id)
  mask=ll.reference_evaluable.eq(True).to_numpy() & ll.label.notna().to_numpy()
  y=ll.loc[mask,'label'].astype(int).to_numpy(); xx=x[mask];ids=idx.loc[mask,'patient_id'].to_numpy()
  elig.append(dict(endpoint=ep,total_patients=len(idx),evaluable=len(y),excluded=len(idx)-len(y),class_counts=json.dumps(pd.Series(y).value_counts().sort_index().to_dict())))
  for site in ['A2','C8','E2','BH']:
   m=joblib.load(ROOT/f'models/{ep}_{site}_model.joblib');assert np.array_equal(m['classes'],np.arange(len(m['classes'])))
   raw=m['raw_model'].predict_proba(xx);cal=calibrated(m,raw)
   for variant,p in [('raw',raw),('calibrated',cal)]:
    assert np.isfinite(p).all() and np.allclose(p.sum(1),1)
    out=dict(endpoint=ep,model_holdout_site=site,variant=variant,n=len(y),**metrics(y,p))
    # Ordinary patient bootstrap preserves observed prevalence variation. Retry absent-class draws.
    rng=np.random.default_rng(20261022);samples=[];attempts=0
    while len(samples)<REPS and attempts<REPS*20:
     attempts+=1;sel=rng.integers(0,len(y),len(y))
     if len(np.unique(y[sel]))!=p.shape[1]:continue
     samples.append(list(metrics(y[sel],p[sel]).values()))
    a=np.array(samples);assert len(a)==REPS
    for j,name in enumerate(metrics(y,p)):
     out[name+'_lo'],out[name+'_hi']=np.quantile(a[:,j],[.025,.975]).tolist()
    out['bootstrap_valid']=len(a);out['bootstrap_attempts']=attempts;rows.append(out)
    for j,pid in enumerate(ids):
     rr=dict(endpoint=ep,model_holdout_site=site,variant=variant,patient_id=pid,label=int(y[j]))
     rr.update({f'prob_class_{c}':float(v) for c,v in enumerate(p[j])});preds.append(rr)
    conf=p[:,1] if p.shape[1]==2 else p.max(1);obs=y if p.shape[1]==2 else (p.argmax(1)==y).astype(int);bi=np.minimum((conf*10).astype(int),9)
    for b in range(10):
     sel=bi==b;bins_out.append(dict(endpoint=ep,model_holdout_site=site,variant=variant,bin=b,n=int(sel.sum()),mean_probability=float(conf[sel].mean()) if sel.any() else None,observed_fraction=float(obs[sel].mean()) if sel.any() else None))
   print(ep,site,'complete',flush=True)
 pd.DataFrame(rows).to_csv(ROOT/'external_metrics_95CI.csv',index=False)
 pd.DataFrame(preds).to_csv(ROOT/'external_predictions.csv',index=False)
 pd.DataFrame(elig).to_csv(ROOT/'endpoint_eligibility.csv',index=False)
 pd.DataFrame(bins_out).to_csv(ROOT/'reliability_bins.csv',index=False)
 summary=dict(status='FROZEN_EXTERNAL_PREDICTIONS_AND_METRICS_COMPLETE_TISSUE_REVIEW_PENDING',patients=len(idx),model_count=20,bootstrap_reps=REPS,bootstrap='ordinary patient percentile bootstrap; absent-class draws retried',training_performed=False,external_calibration_fitted=False,ensemble_created=False,model_selection_performed=False,tumor_region_confirmation=False,same_assayed_block_proven=False,pretraining_independence_established=False,interpretation='All four pre-existing site-holdout models reported separately; no single final all-TCGA model assumed. These overlapping models are not independent replications. Tissue-QC cohort, not tumor-region-confirmed cohort.',multiclass_Brier='unnormalised sum across classes',multiclass_ECE='top-label confidence ECE; binary ECE uses positive probability',label_classes={'binary':{'0':'negative / no qualifying reported mutation','1':'positive / qualifying mutation'},'PAM50':{'0':'LumA','1':'LumB','2':'Basal','3':'Her2'}},inputs_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['inputs','models'] for p in sorted((ROOT/folder).glob('*'))})
 (ROOT/'summary.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':main()
