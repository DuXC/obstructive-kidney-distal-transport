#!/usr/bin/env python3
"""Locked animal-level comparisons; fluorescence and RNA scores are never pooled."""
from pathlib import Path
import json,hashlib,math,sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

P=Path(__file__).resolve().parents[1]
OUT=Path(sys.argv[1]) if len(sys.argv)>1 else P/'05_results'
OUT.mkdir(parents=True,exist_ok=True)
freeze=json.loads((P/'03_protocol/freeze_receipt.json').read_text())
assert hashlib.sha256((P/freeze['file']).read_bytes()).hexdigest()==freeze['sha256']
mods=[m for m in json.loads((P.parent/'03_protocol/modules_v1.json').read_text()) if m['family']=='primary_transport']
gene_names=sorted({g for m in mods for g in m['genes']})
orth=pd.read_csv(P/'02_audit/primary_orthology.tsv',sep='\t')
one={g:d.iloc[0].mouse_symbol for g,d in orth.groupby('human_symbol') if len(d)==1 and d.iloc[0]['type']=='ortholog_one2one'}
meta=pd.read_csv(P/'02_audit/sample_batches.tsv',sep='\t')
meta=meta[meta.condition.ne('Normal')].sort_values('sample').set_index('sample')
meta['scan_day']=meta.scan_date.str.strip().str.split().str[0]
norm=pd.read_csv(P/'01_sources/cy5_normalized_log2.tsv.gz',sep='\t',index_col=0)
features=pd.read_csv(P/'01_sources/target_feature_values.tsv',sep='\t')
features['feature_key']=features.FeatureNum.astype(str)+'|'+features.ProbeName
features=features[features.feature_key.isin(norm.index)&features['sample'].isin(meta.index)].copy()
features['log2_signal']=[norm.loc[k,s] for k,s in zip(features.feature_key,features['sample'])]
flags=['rIsSaturated','rIsFeatNonUnifOL','rIsFeatPopnOL','IsManualFlag']
features['clean']=features[flags].eq(0).all(axis=1)
features['detectable']=features.clean & features.rIsWellAboveBG.eq(1)
features.to_csv(OUT/'target_spot_normalized_and_quality.tsv',sep='\t',index=False)
long=[]
for (sample,mouse),d in features.groupby(['sample','mouse_symbol']):
    for analysis,mask in [('primary',d.clean),('detectable',d.detectable)]:
        probes=d.loc[mask].groupby('ProbeName').log2_signal.median()
        for human in [h for h,m in one.items() if m==mouse]:
            long.append({'sample':sample,'human_symbol':human,'mouse_symbol':mouse,'analysis':analysis,
                'value':probes.median() if len(probes) else np.nan,'usable_probes':len(probes),
                'available_probes':d.ProbeName.nunique(),'clean_spots':int(d.clean.sum()),
                'above_background_spots':int(d.detectable.sum())})
genes=pd.DataFrame(long).sort_values(['analysis','human_symbol','sample'])
genes.to_csv(OUT/'gene_values.tsv',sep='\t',index=False)
scores=[];coverage=[]
for m in mods:
    missing=[g for g in m['genes'] if g not in one]
    coverage.append({'module':m['module'],'members':len(m['genes']),'one_to_one_members':len(m['genes'])-len(missing),'ambiguous_members':';'.join(missing),'testable_mapping':not missing})
    for analysis in ['primary','detectable']:
        wide=genes[genes.analysis.eq(analysis)].pivot(index='sample',columns='human_symbol',values='value').reindex(index=meta.index,columns=m['genes'])
        value=wide.mean(axis=1,skipna=False) if not missing else pd.Series(np.nan,index=meta.index)
        for sample,v in value.items():scores.append({'sample':sample,'module':m['module'],'analysis':analysis,'value':v})
scores=pd.DataFrame(scores)
scores.to_csv(OUT/'program_scores.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(OUT/'program_mapping_coverage.tsv',sep='\t',index=False)

contrasts=[]
for duration in [1,5]:
    for day in [0,10,28]:
        contrasts.append({'id':f'd{duration}_r{day}_vs_sham','duration':duration,'kind':'same_time','day':day,
          'parts':[('RUUO',day,1),('Sham',day,-1)]})
    contrasts.append({'id':f'd{duration}_r28_vs_r0','duration':duration,'kind':'late_vs_acute','day':28,
      'parts':[('RUUO',28,1),('RUUO',0,-1)]})
    contrasts.append({'id':f'd{duration}_sham_adjusted_change','duration':duration,'kind':'difference_in_differences','day':28,
      'parts':[('RUUO',28,1),('Sham',28,-1),('RUUO',0,-1),('Sham',0,1)]})
(OUT/'contrast_definitions.json').write_text(json.dumps(contrasts,indent=2))

def welch(parts):
    ns=[len(a) for a,c in parts]
    if min(ns)<4:return {'status':'insufficient_animals','n_groups':';'.join(map(str,ns))}
    weights=np.array([np.var(a,ddof=1)/len(a)*c*c for a,c in parts])
    effect=float(sum(np.mean(a)*c for a,c in parts))
    variance=weights.sum();df=float(variance**2/sum(v*v/(n-1) for v,n in zip(weights,ns))) if variance>0 else np.nan
    if not np.isfinite(df):return {'status':'zero_or_invalid_variance','n_groups':';'.join(map(str,ns))}
    se=float(np.sqrt(variance));crit=stats.t.ppf(.975,df)
    return {'status':'tested','n_groups':';'.join(map(str,ns)),'effect':effect,'se':se,'df':df,
            'ci_low':effect-crit*se,'ci_high':effect+crit*se,'p':float(2*stats.t.sf(abs(effect/se),df))}

def calculate(values,key,entities,analysis):
    wide=values[values.analysis.eq(analysis)].pivot(index='sample',columns=key,values='value').reindex(meta.index)
    rows=[]
    for entity in entities:
        v=wide[entity] if entity in wide else pd.Series(np.nan,index=meta.index)
        for c in contrasts:
            parts=[]
            for cond,day,coefficient in c['parts']:
                ix=meta.duration_days.eq(c['duration'])&meta.recovery_days.eq(day)&meta.condition.eq(cond)
                parts.append((v[ix].dropna().to_numpy(),coefficient))
            row={'entity':entity,'analysis':analysis,'contrast':c['id'],'duration_days':c['duration'],'recovery_days':c['day'],'kind':c['kind']}
            mapping_ok=entity in one if key=='human_symbol' else all(g in one for g in next(m for m in mods if m['module']==entity)['genes'])
            row.update(welch(parts) if mapping_ok else {'status':'orthology_ambiguous','n_groups':';'.join(str(len(a)) for a,c in parts)})
            row['acquisition_caveat']='era_confounded' if c['kind']=='late_vs_acute' else ('untestable_group_by_era_assumption' if c['kind']=='difference_in_differences' else '')
            rows.append(row)
    d=pd.DataFrame(rows)
    adjusted=multipletests(d.get('p',pd.Series(np.nan,index=d.index)).fillna(1),method='fdr_bh')[1]
    d['fdr']=np.where(d.status.eq('tested'),adjusted,np.nan)
    d['family_slots']=len(d)
    return d

def batch_sensitivity(values,key,entities):
    wide=values[values.analysis.eq('primary')].pivot(index='sample',columns=key,values='value').reindex(meta.index)
    rows=[]
    for entity in entities:
        v=wide[entity] if entity in wide else pd.Series(np.nan,index=meta.index)
        for c in [x for x in contrasts if x['kind']=='same_time']:
            d=meta[meta.duration_days.eq(c['duration'])&meta.recovery_days.eq(c['day'])].copy()
            d['value']=v.reindex(d.index);d=d[d.value.notna()]
            shared=[date for date,t in d.groupby('scan_day') if set(t.condition)=={'RUUO','Sham'}]
            d=d[d.scan_day.isin(shared)]
            n1=int(d.condition.eq('RUUO').sum());n0=int(d.condition.eq('Sham').sum())
            row={'entity':entity,'contrast':c['id'],'duration_days':c['duration'],'recovery_days':c['day'],'n_obstruction':n1,'n_sham':n0,'shared_days':';'.join(sorted(shared))}
            mapping_ok=entity in one if key=='human_symbol' else all(g in one for g in next(m for m in mods if m['module']==entity)['genes'])
            if not mapping_ok:row['status']='orthology_ambiguous'
            elif min(n1,n0)<4:row['status']='insufficient_date_overlap'
            else:
                design=pd.DataFrame({'intercept':1.,'obstruction':d.condition.eq('RUUO').astype(float)},index=d.index)
                design=pd.concat([design,pd.get_dummies(d.scan_day,prefix='date',drop_first=True,dtype=float)],axis=1)
                if len(d)-np.linalg.matrix_rank(design)<4:row['status']='insufficient_residual_df'
                else:
                    fit=sm.OLS(d.value.astype(float),design).fit(cov_type='HC3',use_t=True)
                    ci=fit.conf_int().loc['obstruction']
                    row.update(status='tested',effect=float(fit.params['obstruction']),se=float(fit.bse['obstruction']),df=float(fit.df_resid),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(fit.pvalues['obstruction']))
            rows.append(row)
    df=pd.DataFrame(rows)
    df['fdr']=np.where(df.status.eq('tested'),multipletests(df.get('p',pd.Series(np.nan,index=df.index)).fillna(1),method='fdr_bh')[1],np.nan)
    df['family_slots']=len(df)
    return df

outputs={}
for kind,values,key,entities in [('program',scores,'module',[m['module'] for m in mods]),('gene',genes,'human_symbol',gene_names)]:
    for analysis in ['primary','detectable']:
        d=calculate(values,key,entities,analysis);name=f'{kind}_{analysis}_contrasts'
        d.to_csv(OUT/(name+'.tsv'),sep='\t',index=False);outputs[name]=d
    d=batch_sensitivity(values,key,entities);name=f'{kind}_date_adjusted'
    d.to_csv(OUT/(name+'.tsv'),sep='\t',index=False);outputs[name]=d
summary={name:{'slots':len(d),'tested':int(d.status.eq('tested').sum()),'fdr_below_005':int(d.fdr.lt(.05).sum())} for name,d in outputs.items()}
summary.update(protocol_sha256=freeze['sha256'],animal_count=len(meta),normal_controls_descriptive=4)
(OUT/'results_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print(outputs['program_primary_contrasts'][['entity','contrast','n_groups','effect','ci_low','ci_high','fdr','status']].to_string(index=False))
print(outputs['gene_primary_contrasts'].query('entity in ["SLC12A1","SLC12A3","TRPM6"]')[['entity','contrast','n_groups','effect','ci_low','ci_high','fdr','status']].to_string(index=False))
