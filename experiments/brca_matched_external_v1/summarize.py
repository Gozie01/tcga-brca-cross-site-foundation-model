from pathlib import Path
import sys,json,shutil
import numpy as np,pandas as pd
from run_matched import distributions,REPS
R=Path(__file__).resolve().parent
# Paired raw-vs-calibrated effects, diagnostic only; no probabilities are refitted.
uni=pd.read_csv(R/'inputs/uni_predictions.csv').rename(columns={'model_holdout_site':'site'});hib=pd.read_csv(R/'hibou_external_predictions.csv');rows=[]
for enc,d in [('UNI',uni),('Hibou-B',hib)]:
 for (ep,site),g in d.groupby(['endpoint','site']):
  a=g[g.variant.eq('raw')].sort_values('patient_id');b=g[g.variant.eq('calibrated')].sort_values('patient_id');assert np.array_equal(a.patient_id,b.patient_id);y=a.label.to_numpy(int);k=len(np.unique(y));cols=[f'prob_class_{c}' for c in range(k)];pa=a[cols].to_numpy();pb=b[cols].to_numpy();rng=np.random.default_rng(20261022);draws=[]
  while len(draws)<REPS:
   sel=rng.integers(0,len(y),len(y))
   if len(np.unique(y[sel]))==k:draws.append(np.bincount(sel,minlength=len(y)))
  w=np.vstack([np.ones(len(y)),draws]);aa=distributions(y,pa,w);bb=distributions(y,pb,w)
  for metric in aa:
   delta=bb[metric]-aa[metric];lo,hi=np.percentile(delta[1:],[2.5,97.5]);rows.append(dict(encoder=enc,endpoint=ep,site=site,metric=metric,delta_calibrated_minus_raw=delta[0],ci_lower=lo,ci_upper=hi,valid_bootstrap=REPS))
pd.DataFrame(rows).to_csv(R/'paired_calibration_effects.csv',index=False)
# Original held-out TCGA predictions from the same trained models, for context.
from frozen_metrics import metrics
h=pd.read_csv(R/'hibou_tcga_site_predictions.csv');u=pd.read_csv(R/'inputs/tcga_uni_predictions.csv');context=[]
for (ep,site,var),g in h.groupby(['endpoint','model_holdout_site','variant']):
 g=g.sort_values('patient_id');ref=u[u.endpoint.eq(ep)&u.holdout_site.eq(site)].set_index('patient_id').loc[g.patient_id];y=g.label.to_numpy(int);k=len(np.unique(y));assert np.array_equal(y,ref.label.to_numpy());prefix='prob_raw_class_' if var=='raw' else 'prob_calibrated_class_'
 for enc,p in [('Hibou-B',g[[f'prob_class_{c}' for c in range(k)]].to_numpy()),('UNI',ref[[f'{prefix}{c}' for c in range(k)]].to_numpy())]:context.append(dict(endpoint=ep,site=site,variant=var,encoder=enc,n=len(y),**metrics(y,p)))
pd.DataFrame(context).to_csv(R/'tcga_heldout_context_metrics.csv',index=False)
# Compact summary and figure.
d=pd.read_csv(R/'matched_external_metrics.csv');c=pd.read_csv(R/'paired_encoder_contrasts.csv');a=d[d.metric.eq('auroc')&d.variant.eq('calibrated')]
text='''# Matched external UNI versus Hibou-B comparison

All 20 Hibou-B classifiers were trained using the original TCGA endpoint/site train partitions, with the original training-only C grid, preprocessing, tuning criterion, and seed. Calibration used only the corresponding TCGA validation partition. All models were locked before external evaluation. UNI uses its pre-existing frozen checkpoints and verified predictions. No external tuning, recalibration, model selection, or ensemble was performed.

| Endpoint | External n | UNI calibrated AUROC range | Hibou-B calibrated AUROC range |
|---|---:|---:|---:|
'''
for ep,g in a.groupby('endpoint'):
 u=g[g.encoder.eq('UNI')];h=g[g.encoder.eq('Hibou-B')];text+=f'| {ep} | {int(g.n.iloc[0])} | {u.estimate.min():.3f}–{u.estimate.max():.3f} | {h.estimate.min():.3f}–{h.estimate.max():.3f} |\n'
text+='''
Ranges describe the four source models, not confidence intervals. Per-model intervals and paired differences are provided separately. UNI has higher calibrated AUROC point estimates in all 20 comparisons. Paired Hibou-minus-UNI AUROC intervals exclude zero in the UNI-favouring direction for PAM50 in 4/4 source models, ER in 3/4, PIK3CA in 2/4, and PR in 1/4; none for HER2. These are exploratory, unadjusted intervals across overlapping models and endpoints, not independent confirmatory replications. Both encoders perform weakly for HER2. Probability calibration is not uniformly portable.

## Files

- matched_external_metrics.csv: per-encoder, per-model raw/calibrated metrics and paired-bootstrap-derived marginal intervals.
- paired_encoder_contrasts.csv: Hibou minus UNI paired changes with 95% percentile intervals, 1,000 valid patient resamples. Higher AUROC/AP is better; lower Brier/log loss/ECE is better.
- paired_calibration_effects.csv: calibrated minus raw changes within each encoder/model, using identical patient resamples. No external recalibration is applied.
- calibration_diagnostics.csv: descriptive calibration-in-the-large, joint recalibration intercept and slope. These regression fits diagnose predictions; they are not applied to reported probabilities. PAM50 uses one-vs-rest diagnostics.
- hibou_external_predictions.csv and hibou_tcga_site_predictions.csv: auditable probabilities.
- tcga_heldout_context_metrics.csv: source-domain held-out point estimates for these exact models; not a rerun of all internal five-fold comparisons.
- tcga_fit_audit.csv, models_locked.json, protocol.json: model choices, sizes, source hashes and fixed protocol.
- models and inputs: frozen models and the data needed to reproduce this analysis.

## Outstanding evidence

Both encoders use the same 391 external slides, patient ordering, source-file hashes, patch counts, and coordinate hashes; TCGA Hibou feature hashes and patient indexing pass the supplied freeze. The supplied TCGA Hibou freeze contains no original checkpoint revision or weight hashes. Source code specifies official Hibou processing and float32 inference, and external extraction pins revision 7fdc8cce365ebfd9e7e60ee6ee9a8d4d59388f3c, but identical historical TCGA weights or a cross-cohort numerical replay are not established. Resolve that provenance before a definitive encoder-superiority claim.

Human tumor-content review remains pending. Matching to the molecularly assayed block and encoder pretraining independence are not proven. These computational results therefore remain conditional on those limitations. Missing/equivocal/unmapped endpoint labels are excluded under the frozen reference definitions. No clinical utility or assay-replacement claim is supported.

Binary positive-class ECE is provided alongside the original manuscript top-label ECE (ece_10bin); they are distinct definitions. Multiclass Brier is the unnormalised four-class sum. Patient bootstrap intervals condition on fixed trained models, preserve patient pairing, and do not include training uncertainty.

## Reproduce

Install requirements.txt in a separate environment. Run python run_matched.py, then python summarize.py. No GPU, slide downloads, or feature extraction are required. A clean run fits Hibou models on TCGA only; an existing package reuses the bundled models. Existing UNI results are supplied and were verified against stored TCGA predictions in the earlier external evaluation.
'''
(R/'RESULTS.md').write_text(text)
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,5,figsize=(16,4),sharey=True)
for ax,(ep,g) in zip(axes,a.groupby('endpoint')):
 for j,(enc,col) in enumerate([('UNI','#246589'),('Hibou-B','#c27734')]):
  h=g[g.encoder.eq(enc)].set_index('site').loc[['A2','C8','E2','BH']];v=h.estimate.to_numpy();ax.errorbar(np.arange(4)+(j-.5)*.15,v,yerr=np.array([v-h.ci_lower,h.ci_upper-v]),fmt='o',capsize=2,color=col,label=enc,markersize=4)
 ax.set_title(ep.replace('_binary','').replace('_primary_4class',''));ax.set_xticks(range(4),['A2','C8','E2','BH']);ax.set_ylim(.1,1.02);ax.grid(axis='y',alpha=.2);ax.set_xlabel('TCGA holdout model')
axes[0].set_ylabel('External calibrated AUROC (95% CI)');axes[-1].legend(fontsize=8);fig.suptitle('Matched external comparison — tumor review pending');fig.tight_layout();fig.savefig(R/'matched_external_auroc.png',dpi=180);plt.close(fig)
