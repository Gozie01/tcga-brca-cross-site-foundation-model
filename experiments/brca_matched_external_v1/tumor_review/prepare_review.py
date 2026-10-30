"""Prepare existing contact sheets for outcome-blinded review; no image inference."""
import argparse,csv,hashlib,html,json,zipfile
from pathlib import Path
R=Path(__file__).resolve().parent

def prepare(path,out):
 mapping=list(csv.DictReader((R/'review_manifest.csv').open(newline='',encoding='utf-8-sig')))
 out.mkdir(parents=True,exist_ok=True);(out/'contact_sheets').mkdir(exist_ok=True)
 with zipfile.ZipFile(path) as z:
  summary=json.loads(z.read('summary.json'));assert len(summary['slides'])==391
  byindex={str(x['site_slide_index']):x for x in summary['slides']}
  for row in mapping:
   orig=byindex[row['original_slide_index']]
   assert hashlib.sha256(orig['file_name'].encode()).hexdigest()==row['file_name_sha256']
   data=z.read('preview_'+row['original_slide_index']+'.jpg')
   (out/'contact_sheets'/(row['review_id']+'.jpg')).write_bytes(data)
 fields=['review_id','sampled_tumor_present','review_scope','slide_decision','reviewer_name','reviewer_qualification','review_date','reason_or_notes']
 target=out/'review_decisions.csv'
 if not target.exists():
  with target.open('w',newline='',encoding='utf-8-sig') as f:
   w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(dict(review_id=x['review_id'],sampled_tumor_present='UNREVIEWED',review_scope='',slide_decision='PENDING',reviewer_name='',reviewer_qualification='',review_date='',reason_or_notes='') for x in mapping)
 links='\n'.join('<li><a target="_blank" href="contact_sheets/'+html.escape(x['review_id'])+'.jpg">'+html.escape(x['review_id'])+'</a></li>' for x in mapping)
 (out/'review_index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Blinded histology review</title><style>body{font:18px system-ui;max-width:1050px;margin:40px auto;line-height:1.6}ul{columns:5}a{color:#185a76}</style><h1>Histology contact-sheet review</h1><p>Open each image and record findings in review_decisions.csv. These are sampled tissue patches, not complete slides. A negative or uncertain preview requires whole-slide review before exclusion. No molecular labels or model predictions are supplied.</p><ul>'+links+'</ul></html>',encoding='utf-8')
 print('Prepared',len(mapping),'contact sheets in',out.resolve());print('Open review_index.html and complete review_decisions.csv. Existing decisions were preserved.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--qc-zip',type=Path,required=True);p.add_argument('--out',type=Path,default=R/'reviewer_pack');a=p.parse_args();prepare(a.qc_zip,a.out)
