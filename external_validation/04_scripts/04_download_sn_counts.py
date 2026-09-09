from pathlib import Path
import urllib.request,hashlib,json,time
from datetime import datetime
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[1]
fn='GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Counts_03282022.RDS.gz'
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183277/suppl/'+fn
out=P/'02_data'/fn;temp=out.with_suffix(out.suffix+'.part')
req=urllib.request.Request(url,headers={'User-Agent':'research-secondary-analysis/1.0'})
start=time.time();last=start;h=hashlib.sha256()
with urllib.request.urlopen(req,timeout=90) as response:
 expected=response.headers.get('Content-Length')
 with temp.open('wb') as f:
  while block:=response.read(4*1024*1024):
   f.write(block);h.update(block)
   if time.time()-last>20:print('Downloaded',f.tell(),'of',expected,flush=True);last=time.time()
 if expected:assert temp.stat().st_size==int(expected)
temp.rename(out)
r={'url':url,'file':str(out.relative_to(P)),'bytes':out.stat().st_size,'sha256':h.hexdigest(),'retrieved_at':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'elapsed_seconds':round(time.time()-start,2)}
(P/'08_qa/sn_count_download.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)

