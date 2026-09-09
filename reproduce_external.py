#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,os,subprocess,sys,datetime
import pandas as pd
R=Path(__file__).resolve().parent;P=R/'external_validation'
for name,value in json.loads((P/'INPUT_SHA256.json').read_text()).items():
 assert hashlib.sha256((P/name).read_bytes()).hexdigest()==value,name
for d in ['05_results','06_figures','08_qa']:(P/d).mkdir(exist_ok=True)
rscript=os.environ.get('RSCRIPT','Rscript')
for assay in ['sn','sc']:
 print('Reproducing external',assay,'from verified donor aggregates',flush=True)
 with (P/f'08_qa/{assay}_reproduction.log').open('w') as log:
  subprocess.run([rscript,str(P/'04_scripts/07_external_models.R'),str(P),assay],check=True,stdout=log,stderr=subprocess.STDOUT)
subprocess.run([sys.executable,str(P/'04_scripts/09_verify_statistics.py')],check=True,stdout=subprocess.DEVNULL)
names=['program_comparisons.tsv','donor_scores.tsv','member_coverage.tsv','normalization.tsv','gene_results.tsv.gz','fixed_gene_results.tsv','lodo_summary.tsv']
rows=[]
for assay in ['sn','sc']:
 for suffix in names:
  name=f'{assay}_{suffix}';actual=P/'05_results'/name;expected=P/'expected'/name
  a=pd.read_csv(actual,sep='\t');b=pd.read_csv(expected,sep='\t')
  pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-10,atol=1e-10)
  rows.append({'file':name,'rows':len(a),'numeric_match':True,'byte_identical':actual.read_bytes()==expected.read_bytes()})
subprocess.run([sys.executable,str(P/'04_scripts/10_figures.py')],check=True)
figures=[]
for name in ['Fig4','FigS5','FigS6']:
 a=P/'06_figures'/f'{name}.png';b=P/'expected'/f'{name}.png'
 same=a.read_bytes()==b.read_bytes()
 figures.append({'file':a.name,'byte_identical':same})
 # Byte identity is reported, not required across different font/renderer versions.
result={'status':'PASS','completed_at':datetime.datetime.now().astimezone().isoformat(),'scope':'Reproduction from included verified donor aggregates; full raw source reaggregation is separate.','tables':rows,'figures':figures,'original_discovery_analysis':'Unchanged from v0.4.0.'}
(P/'08_qa/reproduction_receipt.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'status':'PASS','tables':len(rows),'byte_identical_tables':sum(x['byte_identical'] for x in rows),'figures':figures},indent=2))

