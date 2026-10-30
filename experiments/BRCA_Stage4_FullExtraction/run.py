"""All frozen external patches, sequential encoders, per-slide restart points."""
import argparse,csv,gc,importlib.metadata,io,json,os,shutil,time,traceback,zipfile
from pathlib import Path
from datetime import datetime
import numpy as np
from core import sha,identity,atomic_json,read_rows,valid_array,checkpoint,commit,pool,atomic_npy
from backend import load
COORD='d282035f0e34b2961b9f6e0f141fc6bddbf1d884e9e7e6d2a4decc6c1e4dc6fe'

def read_patch(wsi,row):
    from PIL import Image
    size=int(row['read_size_level0'])
    im=wsi.read_region((int(row['x_level0']),int(row['y_level0'])),0,(size,size)).convert('RGB')
    if im.size!=(256,256):im=im.resize((256,256),Image.Resampling.LANCZOS)
    return im

def write_csv(path,rows):
    tmp=path.with_suffix('.csv.tmp')
    with tmp.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    os.replace(tmp,path)

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--qc-report',type=Path);a=p.parse_args()
    home=Path(__file__).resolve().parent;out=home/'features_v1';out.mkdir(exist_ok=True)
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f');started=time.perf_counter()
    report={'status':'RUNNING','models':{},'errors':[],'expected_slides':391,'expected_patients':120,'expected_patches':368985,'predictions_computed':False,'tumor_content_confirmed':False,'original_tcga_weight_identity_confirmed':False,'patient_pooling':'equal mean of float32 slide means; compatibility with original UNI predictor pending'}
    try:
        import torch,openslide
        if not torch.cuda.is_available():raise RuntimeError('CUDA required')
        torch.manual_seed(20261019);torch.backends.cudnn.benchmark=False
        torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
        env={k:importlib.metadata.version(k) for k in ['torch','torchvision','timm','transformers','huggingface-hub','openslide-python','Pillow','numpy','pandas']}
        report['environment']=env;report['gpu']=torch.cuda.get_device_name(0)
        pins=json.loads((home/'model_pins.json').read_text());qc=a.qc_report or a.project_root/'BRCA_ExternalValidation_Stage2_FullQC'/'brca_stage2_full_qc_20261021_041052_426646.zip'
        with zipfile.ZipFile(qc) as z:
            raw_summary=z.read('summary.json')
            if __import__('hashlib').sha256(raw_summary).hexdigest()!='045c07aba44391ae76401904c80eabd06e6be220333df5cfc2ced5e98b950b04':raise ValueError('QC summary changed')
            summary=json.loads(raw_summary);slides=summary['slides']
            if summary['coordinates_sha256']!=COORD or len(slides)!=391 or len(set(s['patient_id'] for s in slides))!=120 or sum(s['final_retained'] for s in slides)!=368985:raise ValueError('Unexpected cohort')
            if len(set(s['file_name'] for s in slides))!=391:raise ValueError('Duplicate slides')
            for s in slides:read_rows(z,s)
            report['coordinates_sha256']=COORD
            if not (out/'run_lock.json').exists() and shutil.disk_usage(out).free<5*1024**3:raise RuntimeError('At least 5 GiB free space required')
            source_code={n:sha(home/n) for n in ['run.py','core.py','backend.py','model_pins.json','pilot_reference.zip']}
            lock={'code':source_code,'coordinates':COORD,'environment':env,'gpu':report['gpu'],'batch_size':8}
            lockpath=out/'run_lock.json'
            if lockpath.exists() and json.loads(lockpath.read_text())!=lock:raise ValueError('Run environment/code changed; existing checkpoints will not be mixed')
            if not lockpath.exists():atomic_json(lockpath,lock)
            registry_path=out/'source_registry.json'
            source_registry=json.loads(registry_path.read_text()) if registry_path.exists() else {}
            with zipfile.ZipFile(home/'pilot_reference.zip') as pilot:
                reference_report=json.loads(pilot.read('report.json'));pilot_rows=list(csv.DictReader(io.StringIO(pilot.read('pilot_index.csv').decode('utf-8'))))
                for name in ['UNI','Hibou']:
                    folder=out/name;folder.mkdir(exist_ok=True);model_start=time.perf_counter()
                    report['models'][name]={'status':'RUNNING','completed_slides':0,'completed_patches':0}
                    try:
                        encode,settings,dim=load(name,pins[name])
                        if identity(settings)!=identity(reference_report['models'][name]['settings']):raise ValueError('Preprocessing differs from pilot')
                        # Repeat the complete published pilot before processing this encoder.
                        probe_images=[]
                        for file in dict.fromkeys(r['file_name'] for r in pilot_rows):
                            s=next(s for s in slides if s['file_name']==file)
                            with openslide.OpenSlide(s['source_path']) as wsi:
                                probe_images.extend(read_patch(wsi,r) for r in pilot_rows if r['file_name']==file)
                        actual=np.concatenate([encode(probe_images[i:i+8]) for i in range(0,len(probe_images),8)])
                        expected=np.load(io.BytesIO(pilot.read(name+'_pilot_features.npy')),allow_pickle=False)
                        delta=float(np.max(np.abs(actual-expected)))
                        if not np.isfinite(actual).all() or not np.array_equal(actual,expected):raise ValueError('Pilot parity failed; max absolute difference '+str(delta))
                        del probe_images,actual,expected
                        model_identity=identity({'lock':lock,'model':pins[name],'settings':settings})
                        records=[];reused=0;torch.cuda.reset_peak_memory_stats()
                        for number,s in enumerate(slides,1):
                            slide_start=time.perf_counter();path=Path(s['source_path']);stat=path.stat()
                            if stat.st_size!=s['source_size']:raise ValueError('Slide size differs from QC: '+s['file_name'])
                            print(name,str(number)+'/391 hashing/checking '+s['file_name'],flush=True)
                            source_hash=sha(path)
                            if s['file_name'] in source_registry and source_registry[s['file_name']]!=source_hash:raise ValueError('Slide bytes changed across encoders/restarts: '+s['file_name'])
                            if s['file_name'] not in source_registry:
                                source_registry[s['file_name']]=source_hash;atomic_json(registry_path,source_registry)
                            rows=read_rows(z,s);token=identity(s['file_name'])[:20]
                            expected_id=identity({'model':model_identity,'source_sha256':source_hash,'coordinate_sha256':s['checkpoint_sha256'],'rows':len(rows)})
                            m=checkpoint(folder,token,expected_id,dim)
                            was_reused=m is not None
                            if m is not None and (m['rows']!=len(rows) or m['file_name']!=s['file_name'] or m['patient_id']!=s['patient_id']):raise ValueError('Checkpoint index identity changed')
                            if m is None:
                                values=[]
                                with openslide.OpenSlide(str(path)) as wsi:
                                    for i in range(0,len(rows),8):
                                        values.append(encode([read_patch(wsi,r) for r in rows[i:i+8]]))
                                arr=np.concatenate(values)
                                if not valid_array(arr,len(rows),dim):raise ValueError('Invalid extracted features')
                                if path.stat().st_mtime_ns!=stat.st_mtime_ns or path.stat().st_size!=stat.st_size:raise ValueError('Slide changed during extraction')
                                m=commit(folder,token,arr,expected_id,{'token':token,'file_name':s['file_name'],'patient_id':s['patient_id'],'source_sha256':source_hash,'coordinate_sha256':s['checkpoint_sha256'],'seconds_including_hash_read_inference':time.perf_counter()-slide_start})
                            else:reused+=1
                            records.append(m)
                            report['models'][name].update(completed_slides=len(records),completed_patches=sum(m['rows'] for m in records),reused_slides=reused)
                            atomic_json(out/'progress.json',report)
                            print(name,str(number)+'/391 complete; patches='+str(m['rows'])+(' (reused)' if was_reused else ''),flush=True)
                        sf,si,pf,patients=pool(folder,records,dim)
                        atomic_npy(folder/'slide_features.npy',sf);atomic_npy(folder/'patient_features_equal_slide.npy',pf)
                        write_csv(folder/'slide_index.csv',si);write_csv(folder/'patient_index.csv',[{'patient_id':x} for x in patients])
                        write_csv(folder/'checkpoint_inventory.csv',[{k:m[k] for k in ['token','file_name','patient_id','rows','source_sha256','coordinate_sha256','feature_sha256','identity']} for m in records])
                        report['models'][name].update(status='COMPLETE',dimension=dim,patient_count=len(patients),pilot_max_absolute_difference=delta,revision=pins[name]['revision'],settings=settings,model_identity=model_identity,elapsed_seconds_this_invocation=time.perf_counter()-model_start,peak_allocated_gpu_bytes=torch.cuda.max_memory_allocated(),inventory_sha256=sha(folder/'checkpoint_inventory.csv'))
                        atomic_json(folder/'summary.json',report['models'][name])
                    finally:
                        if 'encode' in locals():del encode
                        gc.collect();torch.cuda.empty_cache()
        if any(m.get('completed_slides')!=391 or m.get('completed_patches')!=368985 or m.get('status')!='COMPLETE' for m in report['models'].values()) or len(report['models'])!=2:raise ValueError('Incomplete extraction')
        report['status']='FULL_EXTRACTION_COMPLETE_REVIEW_PENDING'
    except Exception as e:
        report['status']='INCOMPLETE';report['errors'].append(type(e).__name__+': '+str(e));traceback.print_exc()
    finally:
        report['elapsed_seconds_this_invocation']=time.perf_counter()-started
        atomic_json(out/'summary.json',report)
        target=home/('brca_full_extraction_report_'+stamp+'.zip')
        with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
            for n in ['summary.json','run_lock.json','progress.json','source_registry.json']:
                if (out/n).is_file():z.write(out/n,n)
            for name in ['UNI','Hibou']:
                folder=out/name
                if not folder.exists():continue
                for f in sorted(folder.glob('*.json')):z.write(f,name+'/'+f.name)
                for n in ['slide_features.npy','patient_features_equal_slide.npy','slide_index.csv','patient_index.csv','checkpoint_inventory.csv']:
                    if (folder/n).is_file():z.write(folder/n,name+'/'+n)
        print('UPLOAD REPORT ZIP:',target,flush=True)
        print('Full patch features remain in:',out,flush=True)
    return 0 if report['status']=='FULL_EXTRACTION_COMPLETE_REVIEW_PENDING' else 1

if __name__=='__main__':
    # The OS releases this lock even after a crash; concurrent writers are refused.
    lockdir=Path(__file__).resolve().parent/'features_v1';lockdir.mkdir(exist_ok=True)
    with (lockdir/'process.lock').open('a+b') as handle:
        handle.seek(0);handle.write(b'0');handle.flush();handle.seek(0)
        if os.name=='nt':
            import msvcrt
            msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        raise SystemExit(main())

