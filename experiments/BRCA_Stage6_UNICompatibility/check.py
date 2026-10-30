"""Offline diagnostic: compare regenerated TCGA embeddings with stored references."""
import argparse,csv,hashlib,io,json,os,time,traceback,zipfile
from pathlib import Path
from datetime import datetime
os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1';os.environ['HF_HUB_DISABLE_TELEMETRY']='1'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()

def relocate(value,base):
    suffix=value.replace('\\','/').split('/TCGA_BRCA_Mutation_Audit/')[-1]
    p=base/'TCGA_BRCA_Mutation_Audit'/suffix
    return p if p.exists() else Path(value)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--wsi-root',type=Path,default=Path('D:/TCGA_BRCA_WSI'));a=ap.parse_args()
    home=Path(__file__).resolve().parent;out=home/('compatibility_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir()
    report={'status':'RUNNING','sources':{},'errors':[],'samples_per_source':32,'network_access_used':False,'training_run':False,'external_predictions_computed':False,'note':'Finite sample agreement supports numerical compatibility; it cannot prove historical checkpoint identity or pretraining independence.'}
    try:
        import numpy as np,torch,openslide,timm
        from PIL import Image
        from timm.data import create_transform,resolve_model_data_config
        if not torch.cuda.is_available():raise RuntimeError('CUDA unavailable')
        torch.manual_seed(20261019);torch.backends.cudnn.benchmark=False;torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
        sources=json.loads((home/'tcga_sources.json').read_text());selected={};wanted=set()
        for name,meta in sources.items():
            index=relocate(meta['index_file'],a.project_root);matrix=relocate(meta['embeddings_file'],a.project_root)
            print('Verifying stored TCGA source:',name,flush=True)
            if sha(index)!=meta['index_sha256'] or sha(matrix)!=meta['embeddings_sha256']:raise ValueError(name+': stored embedding/index checksum mismatch')
            with index.open(encoding='utf-8-sig',newline='') as f:
                rows=[]
                for r in csv.DictReader(f):
                    rows.append(r)
                    if len(rows)==32:break
            if len(rows)!=32:raise ValueError('Fewer than 32 reference patches')
            arr=np.load(matrix,mmap_mode='r',allow_pickle=False)
            if arr.shape!=tuple(meta['embedding_matrix_shape']) or arr.dtype!=np.float32:raise ValueError('Reference matrix shape/dtype changed')
            reference=np.asarray(arr[[int(r['embedding_index']) for r in rows]],dtype=np.float32)
            selected[name]=(meta,rows,reference);wanted.update(r['file_name'] for r in rows)
        found={n:[] for n in wanted}
        for directory,folders,files in os.walk(a.wsi_root,followlinks=False):
            folders[:]=[f for f in folders if not Path(directory,f).is_symlink()]
            for n in files:
                if n in found:found[n].append(Path(directory,n))
        if any(len(v)!=1 for v in found.values()):raise ValueError('Missing or duplicate reference slides: '+str({k:len(v) for k,v in found.items() if len(v)!=1}))
        roots=[Path.home()/'.cache/huggingface/hub']
        if os.environ.get('HF_HOME'):roots.insert(0,Path(os.environ['HF_HOME'])/'hub')
        for key in ['HF_HUB_CACHE','HUGGINGFACE_HUB_CACHE']:
            if os.environ.get(key):roots.insert(0,Path(os.environ[key]))
        revision='b55a5ec6cade1a39edfe6534189a9b8ca7a022f0'
        snapshot=next((r/'models--MahmoodLab--uni'/'snapshots'/revision for r in roots if (r/'models--MahmoodLab--uni'/'snapshots'/revision).is_dir()),None)
        if snapshot is None:raise FileNotFoundError('Pinned UNI snapshot missing')
        pins={'config.json':'c14c59fb26102bdc95a23a67fb42d3dd04d4d5763ef526cbddf467819a05166e','pytorch_model.bin':'56ef09b44a25dc5c7eedc55551b3d47bcd17659a7a33837cf9abc9ec4e2ffb40'}
        for n,h in pins.items():
            if sha(snapshot/n)!=h:raise ValueError('Changed cached UNI file: '+n)
        model=timm.create_model('hf-hub:MahmoodLab/uni@'+revision,pretrained=True,init_values=1e-5,dynamic_img_size=True).to('cuda').eval()
        config=resolve_model_data_config(model);transform=create_transform(**config,is_training=False)
        report.update(revision=revision,verified_cached_files=pins,torch=torch.__version__,timm=timm.__version__,gpu=torch.cuda.get_device_name(0),data_config=config)
        def encode(images,batch_size):
            result=[]
            for i in range(0,len(images),batch_size):
                batch=torch.stack([transform(im) for im in images[i:i+batch_size]]).to('cuda')
                with torch.inference_mode(),torch.autocast(device_type='cuda',dtype=torch.float16):v=model(batch)
                result.append(v.float().cpu().numpy())
            return np.concatenate(result)
        for name,(meta,rows,reference) in selected.items():
            try:
                if json.dumps(config,sort_keys=True)!=json.dumps(meta['preprocessing'],sort_keys=True):raise ValueError('Original preprocessing config differs')
                images=[]
                for r in rows:
                    with openslide.OpenSlide(str(found[r['file_name']][0])) as wsi:
                        size=int(r['read_size_level0']);im=wsi.read_region((int(r['x_level0']),int(r['y_level0'])),0,(size,size)).convert('RGB')
                        if im.size!=(256,256):im=im.resize((256,256),Image.Resampling.LANCZOS)
                        images.append(im)
                original_batch=encode(images,32);external_batch=encode(images,8)
                def compare(values):
                    if not np.isfinite(values).all() or values.shape!=reference.shape:raise ValueError('Invalid regenerated features')
                    cosine=(values.astype(np.float64)*reference).sum(1)/(np.linalg.norm(values.astype(np.float64),axis=1)*np.linalg.norm(reference.astype(np.float64),axis=1))
                    return {'exact_equal':bool(np.array_equal(values,reference)),'max_absolute_difference':float(np.max(np.abs(values-reference))),'relative_l2':float(np.linalg.norm(values-reference)/np.linalg.norm(reference)),'minimum_cosine_similarity':float(cosine.min())}
                report['sources'][name]={'status':'COMPARISON_COMPLETE','reference_shape':list(reference.shape),'original_batch32':compare(original_batch),'external_batch8':compare(external_batch),'file_names':sorted(set(r['file_name'] for r in rows)),'reference_embeddings_sha256':meta['embeddings_sha256'],'reference_index_sha256':meta['index_sha256']}
                np.save(out/(name+'_reference.npy'),reference,allow_pickle=False);np.save(out/(name+'_regenerated_batch32.npy'),original_batch,allow_pickle=False);np.save(out/(name+'_regenerated_batch8.npy'),external_batch,allow_pickle=False)
                with (out/(name+'_patch_index.csv')).open('w',newline='',encoding='utf-8') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
                print(name,json.dumps(report['sources'][name]),flush=True)
            except Exception as e:report['errors'].append(name+': '+str(e));traceback.print_exc()
        report['status']='COMPARISONS_COMPLETE_REVIEW_REQUIRED' if len(report['sources'])==5 and not report['errors'] else 'INCOMPLETE'
    except Exception as e:report['status']='INCOMPLETE';report['errors'].append(type(e).__name__+': '+str(e));traceback.print_exc()
    finally:
        (out/'report.json').write_text(json.dumps(report,indent=2,default=str),encoding='utf-8')
        with zipfile.ZipFile(out.with_suffix('.zip'),'w',zipfile.ZIP_DEFLATED) as z:
            for f in sorted(out.iterdir()):z.write(f,f.name)
            z.write(home/'tcga_classifier_replay.json','tcga_classifier_replay.json')
        print('UPLOAD REPORT ZIP:',out.with_suffix('.zip'),flush=True)
    return 0 if report['status']=='COMPARISONS_COMPLETE_REVIEW_REQUIRED' else 1

if __name__=='__main__':raise SystemExit(main())
