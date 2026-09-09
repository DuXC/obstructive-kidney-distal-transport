from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,concurrent.futures,time,gzip
from datetime import datetime
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[1]
def fetch(url,path):
 start=time.time()
 req=urllib.request.Request(url,headers={'User-Agent':'research-secondary-analysis/1.0'})
 with urllib.request.urlopen(req,timeout=90) as response:
  expected=response.headers.get('Content-Length')
  tmp=path.with_suffix(path.suffix+'.part')
  with tmp.open('wb') as f:
   while block:=response.read(1024*1024): f.write(block)
  if expected:assert tmp.stat().st_size==int(expected)
  tmp.rename(path)
  row={'url':url,'file':str(path.relative_to(P)),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'retrieved_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'elapsed_seconds':round(time.time()-start,2)}
  print(json.dumps(row),flush=True)
  return row
urls=[]
for acc in ['GSE183277','GSE183276','GSE180394']:
 stem=acc[:-3]+'nnn'
 urls.append((f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/soft/{acc}_family.soft.gz',P/'01_sources'/f'{acc}_family.soft.gz'))
for acc,assay in [('GSE183277','snCv3'),('GSE183276','scCv3')]:
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/{acc}/suppl/'
 for suffix in ['Metadata_03282022.txt.gz','Metadata_Field_Descriptions.txt.gz']:
  fn=f'{acc}_Kidney_Healthy-Injury_Cell_Atlas_{assay}_{suffix}'
  urls.append((base+fn,P/'01_sources'/fn))
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 futs=[ex.submit(fetch,url,path) for url,path in urls]
 for fut in concurrent.futures.as_completed(futs):
  try:rows.append(fut.result())
  except Exception as e:rows.append({'error':str(e)});print(str(e),flush=True)
(P/'08_qa/metadata_download_receipt.json').write_text(json.dumps(rows,indent=2))

