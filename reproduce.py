#!/usr/bin/env python3
from pathlib import Path
import argparse,datetime,hashlib,json,os,subprocess,sys
import pandas as pd
from pandas.testing import assert_frame_equal
ROOT=Path(__file__).resolve().parent
PHASE='09_manuscript_development_20260909'
ap=argparse.ArgumentParser();ap.add_argument('--skip-figures',action='store_true');args=ap.parse_args()
for folder in ['08_logs','05_results','02_data/derived',f'{PHASE}/02_analysis',f'{PHASE}/05_figures',f'{PHASE}/07_qa']:(ROOT/folder).mkdir(parents=True,exist_ok=True)
for rel,expected in json.loads((ROOT/'input_sha256.json').read_text()).items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==expected,rel
commands=[([os.environ.get('RSCRIPT','Rscript'),'04_scripts/04_donor_models.R'],'models'),([sys.executable,'04_scripts/05_context_and_verify.py'],'context'),([sys.executable,f'{PHASE}/06_scripts/01_member_influence.py'],'member_influence')]
commands.append(([sys.executable,f'{PHASE}/06_scripts/10_marker_donor_balance.py'],'donor_marker_balance'))
if not args.skip_figures:commands.append(([sys.executable,f'{PHASE}/06_scripts/04_publication_figures.py'],'figures'))
if not args.skip_figures:commands.append(([sys.executable,'10_manuscript_refinement_20260909/06_scripts/04_figures.py'],'publication_v0_3'))
for cmd,label in commands:
    print('Running '+label,flush=True)
    with (ROOT/'08_logs'/('reproduce_'+label+'.log')).open('w') as log:subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
checks=[]
for rel in json.loads((ROOT/'expected/paths.json').read_text()):
    a=ROOT/rel;b=ROOT/'expected'/rel
    actual=pd.read_csv(a,sep='\t');expected=pd.read_csv(b,sep='\t')
    assert_frame_equal(actual,expected,check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-10)
    checks.append(dict(path=rel,rows=len(actual),columns=len(actual.columns),numeric_and_text_match=True,byte_identical=hashlib.sha256(a.read_bytes()).digest()==hashlib.sha256(b.read_bytes()).digest()))
receipt=dict(status='PASS',completed_at=datetime.datetime.now().astimezone().isoformat(),entry='verified donor aggregates',input_hashes_verified=True,tables=checks,figures_regenerated=not args.skip_figures)
(ROOT/'reproduction_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
