"""Matched external UNI/Hibou comparison, fixed TCGA-only training protocol."""
from pathlib import Path
import hashlib,json,warnings,sys
import numpy as np,pandas as pd,joblib
from scipy.special import expit,softmax
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV,StratifiedKFold
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
from frozen_metrics import calibrated,metrics
from bootstrap_metrics import metric_distributions,calibration_diagnostics
R=Path(__file__).resolve().parent;SEED=20260905;REPS=1000

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fitcal(p,y,classes):
 if len(classes)==2:
  q=np.clip(p[:,1],1e-7,1-1e-7);c=LogisticRegression(C=1e6,solver='lbfgs',max_iter=5000,random_state=SEED).fit(np.log(q/(1-q)).reshape(-1,1),y)
  return dict(method='Platt_sigmoid_validation',parameters=dict(coef=c.coef_.tolist(),intercept=c.intercept_.tolist()))
 temps=np.logspace(-1,1,401);logp=np.log(np.clip(p,1e-7,1));loss=[-np.log(np.clip(softmax(logp/t,axis=1)[np.arange(len(y)),y],1e-7,1)).mean() for t in temps]
 return dict(method='temperature_validation',parameters=dict(temperature=float(temps[np.argmin(loss)])))

def distributions(y,p,w):
 d=metric_distributions(y,p,w)
 # Retain original top-label ECE and add positive-class ECE for consistency with external v1.
 if p.shape[1]==2:
  bins=np.minimum((p[:,1]*10).astype(int),9)
  d['positive_ece10']=sum(np.abs(w@((bins==b)*(y-p[:,1]))) for b in range(10))/w.sum(1)
 return d

def main():
 x=np.load(R/'inputs/tcga_hibou.npy');xi=pd.read_csv(R/'inputs/tcga_index.csv');split=pd.read_csv(R/'inputs/site_splits.csv');ex=np.load(R/'inputs/external_hibou.npy');ei=pd.read_csv(R/'inputs/external_index.csv');labs=pd.read_csv(R/'inputs/labels.csv');uni=pd.read_csv(R/'inputs/uni_predictions.csv');old=pd.read_csv(R/'inputs/tcga_uni_predictions.csv')
 assert x.shape==(345,768) and ex.shape==(120,768) and np.isfinite(x).all() and np.isfinite(ex).all()
 assert xi.patient_id.is_unique and ei.patient_id.is_unique and not split.duplicated(['endpoint','holdout_site','patient_id']).any()
 assert np.array_equal(xi.set_index('feature_row').loc[split.feature_row,'patient_id'].to_numpy(),split.patient_id.to_numpy())
 assert not set(xi.patient_id)&set(ei.patient_id)
 protocol={'status':'LOCKED_BEFORE_HIBOU_FITS','comparison':'All four existing site-holdout models, no external model selection or ensemble','training':'TCGA frozen train partition only; StandardScaler + balanced L2 LR; C grid .001,.01,.1,1,10,100; stratified min(4,minority_n) CV; binary AP / multiclass negative log loss; seed 20260905; max_iter 5000','calibration':'TCGA frozen validation only; original binary Platt / PAM50 temperature protocol','external':'Same endpoint-evaluable patients and labels as external UNI v1','bootstrap':'1000 paired ordinary patient resamples; absent-class draws rejected; exploratory percentile CIs without multiplicity adjustment','known_limits':['tumor review pending','historical TCGA Hibou model revision/hash unavailable in supplied freeze; source code uses official processor and float32 inference, but cross-cohort numerical replay not demonstrated','encoder pretraining independence unestablished','image-to-assayed-block correspondence unproven'],'code_sha256':sha(Path(__file__)),'input_hashes':{p.name:sha(p) for p in sorted((R/'inputs').glob('*'))}}
 pp=R/'protocol.json'
 if pp.exists():assert json.loads(pp.read_text())==protocol
 else:pp.write_text(json.dumps(protocol,indent=2))
 # Finish all TCGA-only fits and lock models before reading external endpoint labels for evaluation.
 fitrows=[];sitepred=[]
 for ep in sorted(split.endpoint.unique()):
  for site in ['A2','C8','E2','BH']:
   d=split[split.endpoint.eq(ep)&split.holdout_site.eq(site)];tr=d[d.partition.eq('train')];va=d[d.partition.eq('validation')];te=d[d.partition.eq('test')]
   assert not(set(tr.patient_id)&set(va.patient_id)) and not(set(pd.concat([tr,va]).patient_id)&set(te.patient_id))
   y=tr.label.to_numpy(int);yv=va.label.to_numpy(int);fp=R/f'models/{ep}_{site}_hibou.joblib'
   with warnings.catch_warnings():
    warnings.simplefilter('error',ConvergenceWarning)
    if fp.exists():m=joblib.load(fp)
    else:
     pipe=Pipeline([('scaler',StandardScaler()),('classifier',LogisticRegression(class_weight='balanced',solver='lbfgs',max_iter=5000,random_state=SEED))]);gs=GridSearchCV(pipe,{'classifier__C':[.001,.01,.1,1.,10.,100.]},scoring='average_precision' if len(np.unique(y))==2 else 'neg_log_loss',cv=StratifiedKFold(min(4,int(pd.Series(y).value_counts().min())),shuffle=True,random_state=SEED),n_jobs=1,error_score='raise').fit(x[tr.feature_row],y)
     raw=gs.best_estimator_;classes=raw.classes_;m={'raw_model':raw,'classes':classes,'calibration':fitcal(raw.predict_proba(x[va.feature_row]),yv,classes)};joblib.dump(m,fp)
   fitrows.append(dict(endpoint=ep,site=site,train_n=len(tr),validation_n=len(va),site_test_n=len(te),best_C=m['raw_model'].named_steps['classifier'].C,model_sha256=sha(fp)))
   raw=m['raw_model'].predict_proba(x[te.feature_row]);cal=calibrated(m,raw)
   for variant,p in [('raw',raw),('calibrated',cal)]:
    frame=te[['patient_id','label']].copy();frame['endpoint']=ep;frame['model_holdout_site']=site;frame['variant']=variant
    for c in range(p.shape[1]):frame[f'prob_class_{c}']=p[:,c]
    sitepred.append(frame)
   print('TCGA FIT',ep,site,'complete',flush=True)
 pd.DataFrame(fitrows).to_csv(R/'tcga_fit_audit.csv',index=False)
 (R/'models_locked.json').write_text(json.dumps({p.name:sha(p) for p in sorted((R/'models').glob('*'))},indent=2))
 preds=[];points=[];contrasts=[];diagnostics=[]
 for ep in sorted(split.endpoint.unique()):
  ll=labs[labs.endpoint.eq(ep)].set_index('patient_id').reindex(ei.patient_id);mask=ll.reference_evaluable.eq(True)&ll.label.notna();y=ll.loc[mask,'label'].to_numpy(int);ids=ei.loc[mask.to_numpy(),'patient_id'].to_numpy();xx=ex[mask.to_numpy()];k=len(np.unique(y))
  rng=np.random.default_rng(20261022);draws=[]
  while len(draws)<REPS:
   sel=rng.integers(0,len(y),len(y))
   if len(np.unique(y[sel]))==k:draws.append(np.bincount(sel,minlength=len(y)))
  w=np.vstack([np.ones(len(y)),draws])
  for site in ['A2','C8','E2','BH']:
   m=joblib.load(R/f'models/{ep}_{site}_hibou.joblib');assert np.array_equal(m['classes'],np.arange(k));raw=m['raw_model'].predict_proba(xx)
   for variant,hp in [('raw',raw),('calibrated',calibrated(m,raw))]:
    ref=uni[uni.endpoint.eq(ep)&uni.model_holdout_site.eq(site)&uni.variant.eq(variant)].set_index('patient_id').loc[ids];assert np.array_equal(ref.label.to_numpy(),y);up=ref[[f'prob_class_{c}' for c in range(k)]].to_numpy();a=distributions(y,up,w);b=distributions(y,hp,w)
    for encoder,p,dist in [('UNI',up,a),('Hibou-B',hp,b)]:
     assert np.isfinite(p).all() and np.allclose(p.sum(1),1)
     # Meaningful check: vectorized bootstrap point estimates must match sklearn.
     check=metrics(y,p)
     for source,target in [('auroc','AUROC'),('auprc','AP'),('brier','Brier'),('log_loss','LogLoss')]:assert np.isclose(dist[source][0],check[target],atol=1e-12)
     for name,val in dist.items():
      lo,hi=np.percentile(val[1:],[2.5,97.5]);points.append(dict(endpoint=ep,site=site,variant=variant,encoder=encoder,metric=name,n=len(y),estimate=val[0],ci_lower=lo,ci_upper=hi))
     for row in calibration_diagnostics(y,p):diagnostics.append(dict(endpoint=ep,site=site,variant=variant,encoder=encoder,**row))
    for name in a:
     val=b[name]-a[name];lo,hi=np.percentile(val[1:],[2.5,97.5]);contrasts.append(dict(endpoint=ep,site=site,variant=variant,metric=name,n=len(y),delta_hibou_minus_uni=val[0],ci_lower=lo,ci_upper=hi,valid_bootstrap=REPS))
    frame=pd.DataFrame({'patient_id':ids,'label':y,'endpoint':ep,'site':site,'variant':variant})
    for c in range(k):frame[f'prob_class_{c}']=hp[:,c]
    preds.append(frame)
  print('PAIRED EVALUATION',ep,'complete',flush=True)
 pd.concat(preds).to_csv(R/'hibou_external_predictions.csv',index=False);pd.concat(sitepred).to_csv(R/'hibou_tcga_site_predictions.csv',index=False)
 pd.DataFrame(points).to_csv(R/'matched_external_metrics.csv',index=False);pd.DataFrame(contrasts).to_csv(R/'paired_encoder_contrasts.csv',index=False);pd.DataFrame(diagnostics).to_csv(R/'calibration_diagnostics.csv',index=False)
 (R/'completion.json').write_text(json.dumps({'status':'COMPUTATIONAL_COMPARISON_COMPLETE_PROVENANCE_AND_TUMOR_REVIEW_PENDING','model_fits':20,'external_training':False,'external_calibration':False,'bootstrap_reps':REPS,'exploratory':True,'protocol_sha256':sha(pp)},indent=2))
if __name__=='__main__':
 with threadpool_limits(limits=1):main()
