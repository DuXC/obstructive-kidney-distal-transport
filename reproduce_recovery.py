#!/usr/bin/env python3
"""Reproduce fixed animal-level effects from verified normalized fluorescence."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,os,shutil,tempfile,platform
R=Path(__file__).resolve().parent; E=R/'recovery_extension'
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--figures',action='store_true',help='redraw Fig5/FigS7 from regenerated result tables')
ap.add_argument('--check-normalization',action='store_true',help='independently rerun limma on supplied foreground/background matrices in a temporary directory')
a=ap.parse_args()
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
inputs=json.loads((E/'INPUT_SHA256.json').read_text())
for rel,want in inputs.items():
    if sha(R/rel)!=want:raise RuntimeError('Input checksum differs: '+rel)
out=E/'reproduced_results';out.mkdir(exist_ok=True)
log=subprocess.run([sys.executable,str(E/'04_scripts/05_recovery_statistics.py'),str(out)],check=True,capture_output=True,text=True)
expected=json.loads((E/'EXPECTED_RESULT_SHA256.json').read_text())
checks=[]
for name,want in expected.items():
    reference=E/'05_results'/name
    if sha(reference)!=want:raise RuntimeError('Reference checksum differs: '+name)
    identical=sha(out/name)==want
    if not identical:
        if not name.endswith('.tsv'):raise RuntimeError('A regenerated file differs: '+name)
        import pandas as pd
        observed=pd.read_csv(out/name,sep='\t');baseline=pd.read_csv(reference,sep='\t')
        pd.testing.assert_frame_equal(observed,baseline,check_exact=False,rtol=1e-10,atol=1e-10)
    checks.append({'file':name,'byte_identical':identical,'values_and_text_match':True,'numeric_tolerance':1e-10})

normalization=None
if a.check_normalization:
    import numpy as np,pandas as pd
    with tempfile.TemporaryDirectory(prefix='kidney_recovery_norm_') as td:
        T=Path(td)
        for d in ['01_sources','02_audit','09_qa']:(T/d).mkdir()
        for f in ['cy5_foreground.tsv.gz','cy5_background.tsv.gz']:shutil.copy2(E/'01_sources'/f,T/'01_sources'/f)
        subprocess.run([os.getenv('RSCRIPT','Rscript'),str(E/'04_scripts/04_normalize_cy5.R'),str(T)],check=True)
        old=pd.read_csv(E/'01_sources/cy5_normalized_log2.tsv.gz',sep='\t',index_col=0)
        new=pd.read_csv(T/'01_sources/cy5_normalized_log2.tsv.gz',sep='\t',index_col=0)
        assert old.index.equals(new.index) and old.columns.equals(new.columns)
        delta=float(np.max(np.abs(new.to_numpy()-old.to_numpy())))
        if delta>1e-10:raise RuntimeError('Normalization mismatch: '+str(delta))
        normalization={'status':'PASS','shape':list(new.shape),'max_absolute_difference':delta,'tolerance':1e-10,'scope':'Included raw foreground/background matrices; no raw source reacquisition'}
if a.figures:
    subprocess.run([sys.executable,str(E/'04_scripts/06_plot_recovery.py'),str(out)],check=True)
receipt={'status':'PASS_RECOVERY_RELEASE_REPRODUCTION','scope':'Verified normalized fluorescence and source quality/metadata inputs; source raw reacquisition is separate',
         'input_hashes_checked':len(inputs),'python':platform.python_version(),'results':checks,'normalization_check':normalization,'figures_regenerated':a.figures}
v=E/'validation';v.mkdir(exist_ok=True)
(v/'reproduction_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
