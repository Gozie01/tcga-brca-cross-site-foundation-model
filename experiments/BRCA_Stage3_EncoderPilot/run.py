"""Offline diagnostic pilot. Outputs are not full-slide or patient features."""
import argparse, csv, gc, hashlib, io, json, os, time, traceback, zipfile
from pathlib import Path
from datetime import datetime
os.environ['HF_HUB_OFFLINE']='1'
os.environ['TRANSFORMERS_OFFLINE']='1'
os.environ['HF_HUB_DISABLE_TELEMETRY']='1'
COORD_HASH='d282035f0e34b2961b9f6e0f141fc6bddbf1d884e9e7e6d2a4decc6c1e4dc6fe'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()

def select(summary):
    chosen=[]
    for mpp in sorted(set(s['source_mpp_x'] for s in summary['slides'])):
        group=[s for s in summary['slides'] if s['source_mpp_x']==mpp]
        first=min(group,key=lambda s:s['file_name'])
        low=min((s for s in group if s['file_name']!=first['file_name']),key=lambda s:(s['final_retained'],s['file_name']))
        chosen.extend([first,low])
    return chosen

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--qc-report',type=Path);a=p.parse_args()
    out=Path(__file__).resolve().parent/('pilot_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir()
    report={'status':'STARTED','errors':[],'models':{},'full_cohort_extraction':False,'predictions_computed':False,'original_tcga_weight_identity_confirmed':False,'patches_per_slide':32}
    try:
        import numpy as np
        import torch, openslide
        from PIL import Image
        if not torch.cuda.is_available():raise RuntimeError('CUDA unavailable; no CPU substitution permitted')
        torch.manual_seed(20261019);torch.backends.cudnn.benchmark=False
        torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
        report['gpu']=torch.cuda.get_device_name(0);report['torch']=torch.__version__
        qc=a.qc_report or a.project_root/'BRCA_ExternalValidation_Stage2_FullQC'/'brca_stage2_full_qc_20261021_041052_426646.zip'
        images=[];index=[]
        with zipfile.ZipFile(qc) as z:
            summary=json.loads(z.read('summary.json'))
            if summary['coordinates_sha256']!=COORD_HASH:raise ValueError('Unexpected cohort coordinate identity')
            report['qc_coordinates_sha256']=COORD_HASH
            for slide in select(summary):
                data=z.read('slide_checkpoints/'+slide['checkpoint_csv'].replace('\\','/').split('/')[-1])
                if hashlib.sha256(data).hexdigest()!=slide['checkpoint_sha256']:raise ValueError('Coordinate checksum mismatch')
                rows=[r for r in csv.DictReader(io.StringIO(data.decode('utf-8-sig'))) if r['retain_after_new_locked_qc'].lower() in ['true','1','1.0']]
                if len(rows)!=slide['final_retained']:raise ValueError('Retained count mismatch')
                rows=[rows[i] for i in np.linspace(0,len(rows)-1,min(32,len(rows)),dtype=int)]
                print('Reading pilot:',slide['file_name'],len(rows),'patches',flush=True)
                with openslide.OpenSlide(slide['source_path']) as wsi:
                    for r in rows:
                        size=int(r['read_size_level0']);im=wsi.read_region((int(r['x_level0']),int(r['y_level0'])),0,(size,size)).convert('RGB')
                        if im.size!=(256,256):im=im.resize((256,256),Image.Resampling.LANCZOS)
                        images.append(im);index.append({k:r[k] for k in ['patient_id','file_name','patch_id','x_level0','y_level0','read_size_level0','source_mpp_x']})
        with (out/'pilot_index.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(index[0]));w.writeheader();w.writerows(index)
        report['patch_count']=len(images);report['slide_count']=len(set(r['file_name'] for r in index))
        pins=json.loads((Path(__file__).parent/'model_pins.json').read_text())
        for name,pin in pins.items():
            try:
                print('Loading offline',name,flush=True)
                roots=[Path.home()/'.cache'/'huggingface'/'hub']
                if os.environ.get('HF_HOME'):roots.insert(0,Path(os.environ['HF_HOME'])/'hub')
                for key in ['HF_HUB_CACHE','HUGGINGFACE_HUB_CACHE']:
                    if os.environ.get(key):roots.insert(0,Path(os.environ[key]))
                snapshot=next((r/pin['repository']/'snapshots'/pin['revision'] for r in roots if (r/pin['repository']/'snapshots'/pin['revision']).is_dir()),None)
                if snapshot is None:raise FileNotFoundError('Pinned cached snapshot absent')
                for item in pin['files']:
                    if sha(snapshot/item['name'])!=item['sha256']:raise ValueError('Changed cached file: '+item['name'])
                if name=='UNI':
                    import timm
                    from timm.data import create_transform,resolve_model_data_config
                    model=timm.create_model('hf-hub:MahmoodLab/uni@'+pin['revision'],pretrained=True,init_values=1e-5,dynamic_img_size=True).to('cuda').eval()
                    config=resolve_model_data_config(model);transform=create_transform(**config,is_training=False)
                    def encode(ims):
                        batch=torch.stack([transform(im) for im in ims]).to('cuda')
                        with torch.inference_mode(),torch.autocast(device_type='cuda',dtype=torch.float16):v=model(batch)
                        return v.float().cpu().numpy()
                    settings={'data_config':config,'transform':str(transform),'precision':'float16_autocast'};dim=1024
                else:
                    from transformers import AutoImageProcessor,AutoModel
                    processor=AutoImageProcessor.from_pretrained(str(snapshot),trust_remote_code=True,local_files_only=True)
                    model=AutoModel.from_pretrained(str(snapshot),trust_remote_code=True,local_files_only=True,use_safetensors=True).to('cuda').eval()
                    def encode(ims):
                        batch=processor(images=ims,return_tensors='pt').to('cuda')
                        with torch.inference_mode():v=model(**batch).pooler_output
                        return v.float().cpu().numpy()
                    settings={'processor':processor.to_dict(),'precision':'float32_no_autocast'};dim=768
                torch.cuda.reset_peak_memory_stats();torch.cuda.synchronize();start=time.perf_counter();parts=[]
                for i in range(0,len(images),8):
                    parts.append(encode(images[i:i+8]));print(name,min(i+8,len(images)),'/',len(images),flush=True)
                torch.cuda.synchronize();seconds=time.perf_counter()-start;features=np.concatenate(parts)
                if features.shape!=(len(images),dim) or not np.isfinite(features).all() or np.any(np.linalg.norm(features,axis=1)==0):raise ValueError('Invalid feature dimensions or values')
                repeat=encode(images[:8]);delta=float(np.max(np.abs(repeat-features[:8])))
                np.save(out/(name+'_pilot_features.npy'),features,allow_pickle=False)
                report['models'][name]={'status':'PILOT_COMPLETE','revision':pin['revision'],'verified_files':pin['files'],'shape':list(features.shape),'extraction_seconds_excluding_slide_reads':seconds,'patches_per_second':len(images)/seconds,'peak_allocated_gpu_bytes':torch.cuda.max_memory_allocated(),'repeat_batch_max_absolute_difference':delta,'settings':settings,'batch_size':8,'feature_sha256':sha(out/(name+'_pilot_features.npy'))}
            except Exception as e:
                report['models'][name]={'status':'FAILED','error':type(e).__name__+': '+str(e)};report['errors'].append(name+': '+str(e));traceback.print_exc()
            finally:
                if 'encode' in locals():del encode
                if 'model' in locals():del model
                gc.collect();torch.cuda.empty_cache()
        report['status']='PILOT_COMPLETE_REVIEW_PENDING' if not report['errors'] else 'PILOT_INCOMPLETE'
    except Exception as e:report['status']='PILOT_INCOMPLETE';report['errors'].append(type(e).__name__+': '+str(e));traceback.print_exc()
    finally:
        (out/'report.json').write_text(json.dumps(report,indent=2,default=str),encoding='utf-8')
        with zipfile.ZipFile(out.with_suffix('.zip'),'w',zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):z.write(f,f.name)
        print('UPLOAD REPORT ZIP:',out.with_suffix('.zip'),flush=True)
    return 0 if report['status']=='PILOT_COMPLETE_REVIEW_PENDING' else 1

if __name__=='__main__':raise SystemExit(main())
