#!/usr/bin/env python3
"""Post hoc member influence diagnostic on the unchanged primary normalized scores."""
from pathlib import Path
import datetime, hashlib, json
import numpy as np
import pandas as pd
from scipy import stats

DEV=Path(__file__).resolve().parents[1]
ROOT=DEV.parent
freeze=json.loads((ROOT/'03_protocol/freeze_receipt.json').read_text())
assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in freeze['files'].items())
mods=json.loads((ROOT/'03_protocol/modules_v1.json').read_text())
rows=[];contrib=[];summ=[]
def contrast(v,group):
    a=v[group=='UUO'];b=v[group=='Control'];delta=a.mean()-b.mean()
    vv=a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)
    df=vv**2/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1))
    half=stats.t.ppf(.975,df)*np.sqrt(vv)
    return float(delta),float(delta-half),float(delta+half)
for mod in mods:
    if mod['family']=='context':continue
    ct=mod['celltype']
    data=pd.read_csv(ROOT/f'02_data/derived/logcpm_{ct}.tsv.gz',sep='\t',index_col=0)
    group=np.array(['UUO' if x.startswith('UUO') else 'Control' for x in data.columns])
    n0=int(sum(group=='Control'));n1=int(sum(group=='UUO'))
    members=mod['genes'];values=data.loc[members].to_numpy();main=contrast(values.mean(axis=0),group)
    deltas=values[:,group=='UUO'].mean(axis=1)-values[:,group=='Control'].mean(axis=1)
    assert np.isclose(deltas.mean(),main[0])
    for g,d in zip(members,deltas):contrib.append(dict(module=mod['module'],celltype=ct,gene=g,mean_log2cpm_difference=float(d),contribution_to_module_difference=float(d/len(members)),absolute_share=float(abs(d)/np.abs(deltas).sum())))
    if min(n0,n1)<4 or len(members)<3:
        summ.append(dict(module=mod['module'],celltype=ct,n_control=n0,n_uuo=n1,n_members=len(members),main_effect=main[0],n_omissions=0,status='TWO_MEMBER_SET' if len(members)<3 else 'INSUFFICIENT_GROUP_COVERAGE'))
        continue
    effects=[]
    for i,g in enumerate(members):
        d,lo,hi=contrast(np.delete(values,i,axis=0).mean(axis=0),group);effects.append(d)
        rows.append(dict(module=mod['module'],celltype=ct,omitted_gene=g,n_control=n0,n_uuo=n1,n_remaining=len(members)-1,effect=d,ci_low=lo,ci_high=hi,same_direction=bool(np.sign(d)==np.sign(main[0])),status='POST_HOC_DIAGNOSTIC'))
    summ.append(dict(module=mod['module'],celltype=ct,n_control=n0,n_uuo=n1,n_members=len(members),main_effect=main[0],n_omissions=len(effects),n_same_direction=int(sum(np.sign(effects)==np.sign(main[0]))),effect_min=min(effects),effect_max=max(effects),largest_abs_member=members[int(np.argmax(np.abs(deltas)))],largest_abs_share=float(np.max(np.abs(deltas))/np.abs(deltas).sum()),status='POST_HOC_DIAGNOSTIC'))
for name,data in [('member_omission_estimates',rows),('member_contributions',contrib),('member_influence_summary',summ)]:pd.DataFrame(data).to_csv(DEV/f'02_analysis/{name}.tsv',sep='\t',index=False)
(DEV/'07_qa/member_influence_receipt.json').write_text(json.dumps({'status':'PASS','completed_at':datetime.datetime.now().astimezone().isoformat(),'frozen_files_unchanged':True,'source':'existing primary TMM log2 CPM; fixed eligible donors','inference_scope':'post hoc member influence diagnostics, no replacement of primary programs or primary multiplicity'},indent=2))
print(pd.DataFrame(summ).to_string(index=False))
