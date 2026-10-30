import os, hashlib
from pathlib import Path
os.environ['HF_HUB_OFFLINE']='1'
os.environ['TRANSFORMERS_OFFLINE']='1'
os.environ['HF_HUB_DISABLE_TELEMETRY']='1'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()

def load(name,pin):
    import torch
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
    return encode, settings, dim
