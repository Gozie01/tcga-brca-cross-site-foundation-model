from __future__ import annotations
import argparse, hashlib, importlib.util, json, os, platform, time
from pathlib import Path
import joblib, numpy as np, pandas as pd, sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

VERSION="stage10c_site_holdout_domain_shift_v1"; SEED=20260905
C_GRID=[.001,.01,.1,1.,10.,100.]; SITES=("A2","C8","E2","BH")
EXPECTED={"ER_binary":337,"PR_binary":337,"HER2_binary":314,"PIK3CA_binary":345,"PAM50_primary_4class":335}

def load_common(root):
 p=root/"Stage10B_Internal_Prediction_Calibration_v1"/"stage10b_run_internal_prediction_calibration_v1.py"
 if not p.is_file(): raise FileNotFoundError(f"Missing Stage 10B support script: {p}")
 spec=importlib.util.spec_from_file_location("stage10b_common",p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""): h.update(b)
 return h.hexdigest()
def atomic_csv(p,d):
 t=p.with_suffix(p.suffix+".tmp"); d.to_csv(t,index=False,encoding="utf-8-sig"); os.replace(t,p)
def atomic_json(p,o):
 t=p.with_suffix(p.suffix+".tmp"); t.write_text(json.dumps(o,indent=2),encoding="utf-8"); os.replace(t,p)
def softmax(z):
 z=z-z.max(1,keepdims=True); e=np.exp(z); return e/e.sum(1,keepdims=True)
def calibrate(raw_val,y_val,raw_test,classes):
 eps=1e-7
 if len(classes)==2:
  pos=int(np.where(classes==1)[0][0]); pv=np.clip(raw_val[:,pos],eps,1-eps); pt=np.clip(raw_test[:,pos],eps,1-eps)
  lv=np.log(pv/(1-pv)).reshape(-1,1); lt=np.log(pt/(1-pt)).reshape(-1,1)
  c=LogisticRegression(C=1e6,solver="lbfgs",max_iter=5000,random_state=SEED).fit(lv,y_val)
  p=c.predict_proba(lt); return p,{"method":"Platt_sigmoid_validation","parameters":{"coef":c.coef_.tolist(),"intercept":c.intercept_.tolist()}}
 logv=np.log(np.clip(raw_val,eps,1)); logt=np.log(np.clip(raw_test,eps,1)); temps=np.logspace(-1,1,401)
 losses=[]
 for temp in temps:
  p=softmax(logv/temp); losses.append(-np.mean(np.log(np.clip(p[np.arange(len(y_val)),y_val],eps,1))))
 temp=float(temps[int(np.argmin(losses))]); return softmax(logt/temp),{"method":"temperature_validation","parameters":{"temperature":temp}}
def model(seed):
 return Pipeline([("scaler",StandardScaler()),("classifier",LogisticRegression(class_weight="balanced",solver="lbfgs",max_iter=5000,random_state=seed))])
def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--root",type=Path,required=True); ap.add_argument("--bootstrap-reps",type=int,default=1000); ap.add_argument("--n-jobs",type=int,default=-1); a=ap.parse_args(); root=a.root; common=load_common(root)
 s10a=root/"stage10a_cohort_split_freeze_v1"; s10b=root/"stage10b_internal_prediction_calibration_v1"
 fp=s10a/"stage10a_uni_patient_features_v1.npy"; sp=s10a/"stage10a_site_holdout_manifest_v1.csv"; mp=s10a/"stage10a_cohort_split_freeze_metadata_v1.json"; ip=s10b/"stage10b_internal_oof_predictions_v1.csv"
 for p in (fp,sp,mp,ip):
  if not p.is_file(): raise FileNotFoundError(p)
 meta=json.loads(mp.read_text(encoding="utf-8"));
 if meta.get("status")!="LOCKED_COMPLETE": raise RuntimeError("Stage 10A is not locked")
 for p in (fp,sp):
  if sha(p)!=meta["outputs"][p.name]["sha256"]: raise RuntimeError(f"Hash mismatch: {p.name}")
 X=np.load(fp,mmap_mode="r"); splits=pd.read_csv(sp); internal=pd.read_csv(ip)
 out=root/VERSION; ck=out/"split_checkpoints"; ck.mkdir(parents=True,exist_ok=True)
 preds=[]; metrics=[]; boot=[]
 print("="*110); print("STAGE 10C TCGA SITE-HOLDOUT DOMAIN SHIFT"); print("="*110)
 for ep,n in EXPECTED.items():
  for site in SITES:
   d=splits[(splits.endpoint==ep)&(splits.holdout_site==site)].copy(); tr=d[d.partition=="train"]; va=d[d.partition=="validation"]; te=d[d.partition=="test"]
   stem=f"{ep}_{site}"; pp=ck/f"{stem}_predictions.csv"; xp=ck/f"{stem}_metrics.csv"; jp=ck/f"{stem}_model.joblib"; qp=ck/f"{stem}_metadata.json"
   fingerprint=hashlib.sha256(f"{sha(fp)}|{sha(sp)}|{VERSION}|{ep}|{site}".encode()).hexdigest()
   valid=False
   if all(p.is_file() for p in (pp,xp,jp,qp)):
    try:
     q=json.loads(qp.read_text()); valid=q.get("fingerprint")==fingerprint and sha(pp)==q["predictions_sha256"] and sha(xp)==q["metrics_sha256"] and sha(jp)==q["model_sha256"]
    except Exception: valid=False
   if valid: pr=pd.read_csv(pp); mr=pd.read_csv(xp); status="VERIFIED SKIP"
   else:
    rows=lambda z:z.feature_row.astype(int).to_numpy(); ytr=tr.label.astype(int).to_numpy(); yv=va.label.astype(int).to_numpy(); yt=te.label.astype(int).to_numpy(); classes=np.sort(d.label.unique().astype(int)); cv=min(4,int(pd.Series(ytr).value_counts().min()))
    scoring="average_precision" if len(classes)==2 else "neg_log_loss"; started=time.perf_counter()
    gs=GridSearchCV(model(SEED),{"classifier__C":C_GRID},scoring=scoring,cv=StratifiedKFold(cv,shuffle=True,random_state=SEED),n_jobs=a.n_jobs,error_score="raise").fit(X[rows(tr)],ytr)
    raw=gs.best_estimator_; rv=raw.predict_proba(X[rows(va)]); rt=raw.predict_proba(X[rows(te)]); cal,calmeta=calibrate(rv,yv,rt,classes)
    pr=te[["endpoint","holdout_site","patient_id","feature_row","source","tss","label","label_name"]].reset_index(drop=True)
    for i,c in enumerate(classes): pr[f"prob_raw_class_{c}"]=rt[:,i]; pr[f"prob_calibrated_class_{c}"]=cal[:,i]
    rr=common.calculate_metrics(yt,rt,classes); cr=common.calculate_metrics(yt,cal,classes); mr=pd.DataFrame([{"endpoint":ep,"holdout_site":site,"calibration":"raw","train_n":len(tr),"validation_n":len(va),"test_n":len(te),"best_C":gs.best_params_["classifier__C"],"elapsed_seconds":time.perf_counter()-started,**rr},{"endpoint":ep,"holdout_site":site,"calibration":calmeta["method"],"train_n":len(tr),"validation_n":len(va),"test_n":len(te),"best_C":gs.best_params_["classifier__C"],"elapsed_seconds":time.perf_counter()-started,**cr}])
    atomic_csv(pp,pr); atomic_csv(xp,mr); tmp=jp.with_suffix(".tmp"); joblib.dump({"raw_model":raw,"calibration":calmeta,"classes":classes},tmp); os.replace(tmp,jp); atomic_json(qp,{"status":"CHECKPOINT_COMPLETE","fingerprint":fingerprint,"predictions_sha256":sha(pp),"metrics_sha256":sha(xp),"model_sha256":sha(jp)}); status="CHECKPOINT PASS"
   preds.append(pr); metrics.append(mr)
   for _,row in mr.iterrows(): print(f"{ep:<24} {site} {row.calibration:<28} AUROC={row.auroc:.4f} ECE={row.ece_10bin:.4f}")
   print(f"  {status}")
 predictions=pd.concat(preds,ignore_index=True); fold_metrics=pd.concat(metrics,ignore_index=True)
 # Site-specific bootstrap confidence intervals for every endpoint and calibration state.
 for ep in EXPECTED:
  for site in SITES:
   subset=predictions[(predictions.endpoint==ep)&(predictions.holdout_site==site)]
   b=common.bootstrap_intervals(subset,ep,a.bootstrap_reps)
   b.insert(1,"holdout_site",site)
   b.loc[b.calibration.eq("sigmoid_cv"),"calibration"]="temperature_validation" if ep=="PAM50_primary_4class" else "Platt_sigmoid_validation"
   boot.append(b)
 bootstrap=pd.concat(boot,ignore_index=True)
 # Aggregate internal reference and TCGA site-held-out domain gap.
 gap=[]
 for ep in EXPECTED:
  dint=internal[internal.endpoint==ep]; classes=np.sort(dint.label.unique().astype(int))
  for state,prefix in (("raw","prob_raw_class_"),("calibrated","prob_calibrated_class_")):
   pi=np.column_stack([dint[f"{prefix}{c}"].to_numpy(float) for c in classes]); im=common.calculate_metrics(dint.label.to_numpy(int),pi,classes)
   for site in SITES:
    rows=fold_metrics[(fold_metrics.endpoint==ep)&(fold_metrics.holdout_site==site)]; row=rows.iloc[0 if state=="raw" else 1]
    gap.append({"endpoint":ep,"holdout_site":site,"calibration_state":state,"internal_auroc":im["auroc"],"external_auroc":row.auroc,"auroc_gap_external_minus_internal":row.auroc-im["auroc"],"internal_auprc":im["auprc"],"external_auprc":row.auprc,"auprc_gap_external_minus_internal":row.auprc-im["auprc"]})
 domain=pd.DataFrame(gap)
 po=out/"stage10c_site_holdout_predictions_v1.csv"; mo=out/"stage10c_site_holdout_metrics_v1.csv"; bo=out/"stage10c_site_holdout_bootstrap_intervals_v1.csv"; go=out/"stage10c_domain_gap_v1.csv"; jo=out/"stage10c_site_holdout_metadata_v1.json"
 atomic_csv(po,predictions); atomic_csv(mo,fold_metrics); atomic_csv(bo,bootstrap); atomic_csv(go,domain)
 expected_test_rows=int(splits[(splits.endpoint.isin(EXPECTED))&(splits.holdout_site.isin(SITES))&(splits.partition=="test")].shape[0])
 gates={"twenty_checkpoints":len(list(ck.glob("*_metadata.json")))==20,"two_states_per_experiment":fold_metrics.groupby(["endpoint","holdout_site"]).size().eq(2).all(),"test_predictions_complete":len(predictions)==expected_test_rows,"no_duplicate_test_predictions":not predictions.duplicated(["endpoint","holdout_site","patient_id"]).any(),"bootstrap_intervals_present":not bootstrap.empty and bootstrap.groupby(["endpoint","holdout_site","calibration"]).ngroups==40,"metrics_finite":np.isfinite(fold_metrics[["auroc","auprc","brier","log_loss","ece_10bin"]].to_numpy()).all()}
 if not all(bool(v) for v in gates.values()): raise RuntimeError(gates)
 atomic_json(jo,{"pipeline_version":VERSION,"status":"LOCKED_COMPLETE","holdout_sites":list(SITES),"endpoint_counts":EXPECTED,"expected_test_rows":expected_test_rows,"bootstrap_reps":a.bootstrap_reps,"calibration":{"binary":"Platt sigmoid on frozen validation","multiclass":"temperature scaling on frozen validation"},"software":{"python":platform.python_version(),"sklearn":sklearn.__version__},"validation":{k:bool(v) for k,v in gates.items()},"outputs":{p.name:{"path":str(p),"sha256":sha(p)} for p in (po,mo,bo,go)}})
 print("="*110); print("STAGE 10C SITE-HOLDOUT AUDIT"); print("Status                         : LOCKED_COMPLETE"); print("Experiments                    : 20"); [print(f"{k:<47}: PASS") for k in gates]; print(f"Metadata                       : {jo}"); print("="*110)
if __name__=="__main__": main()
