"""Retrieve public source metadata into a separate directory for a new audit."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile

URLS={
 'cohort_builder.csv':'https://pathdb.cancerimagingarchive.net/system/files/collectionmetadata/202401/cohort_builder_v1_01-16-2024.csv',
 'archived_specimens.zip':'https://raw.githubusercontent.com/ImagingDataCommons/idc-wsi-conversion/main/MetadataFromTCIAJSONAPIDownload.zip',
 'cbio_gene_matrix_resolved.txt':'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_cptac_2020/data_gene_panel_matrix.txt'
}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);provenance=[]
    for name,url in URLS.items():
        target=a.out/name
        request=urllib.request.Request(url,headers={'User-Agent':'BRCA-CPTAC-Research-Audit/1.0'})
        with urllib.request.urlopen(request,timeout=60) as response,target.open('wb') as f:
            while block:=response.read(1024*1024):f.write(block)
        provenance.append({'url':url,'file':name,'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
        print(name,target.stat().st_size,flush=True)
    with (a.out/'cohort_builder.csv').open(encoding='utf-8',errors='replace',newline='') as f:
        reader=csv.DictReader(f);columns=reader.fieldnames;rows=[r for r in reader if r['collection']=='CPTAC-BRCA']
    with (a.out/'pathdb_cptac_brca_inventory.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    with zipfile.ZipFile(a.out/'archived_specimens.zip') as z:
        (a.out/'archived_BR.json').write_bytes(z.read('MetadataFromTCIAJSONAPIDownload_CPTAC2/BR.json'))
    (a.out/'source_provenance.json').write_text(json.dumps(provenance,indent=2))
    print('New metadata snapshot retrieved. Do not replace the existing acquisition/label lock without a new version.')
