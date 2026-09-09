"""Quantitative manuscript figures from the unchanged donor-level results."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image

DEV=Path(__file__).resolve().parents[1];ROOT=DEV.parent;R=ROOT/'05_results';OUT=DEV/'05_figures'
OUT.mkdir(parents=True,exist_ok=True);(DEV/'00_admin').mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'lines.linewidth':1,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False,'savefig.facecolor':'white'})
COL={'Control':'#446D87','UUO':'#A94F29'}
DON=[f'Control{x}' for x in range(1,8)]+[f'UUO{x}' for x in range(1,6)]
CT=['TAL','DCT','PC','IC-A','IC-B'];PRIMARY=['TAL_salt','TAL_CaMg','DCT_NaCl','DCT_Mg','DCT_Ca']
SHORT={'TAL_salt':'TAL salt','TAL_CaMg':'TAL Ca/Mg','DCT_NaCl':'DCT NaCl','DCT_Mg':'DCT Mg','DCT_Ca':'DCT Ca','PC_water':'PC water','PC_NaK':'PC Na/K','ICA_acid':'IC-A acid','ICB_base':'IC-B base'}
cmp=pd.read_csv(R/'module_comparisons_primary_run.tsv',sep='\t')
scores=pd.read_csv(R/'module_scores_all_scenarios.tsv',sep='\t');main=scores[scores.scenario.eq('main')]
data_manifest=[]
def source(name,files,claim):data_manifest.append(dict(figure=name,source_files=files,claim=claim,backend='matplotlib',archetype='quantitative grid'))
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=12)
    ax.text(-.13,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=11,va='bottom')
def leg(ax,loc='best',ncol=1):
    ax.legend(handles=[Line2D([],[],marker='o' if g=='Control' else '^',linestyle='',color=COL[g],label='RCC reference' if g=='Control' else 'UTUC + obstruction',markersize=4) for g in COL],loc=loc,fontsize=8,ncol=ncol)
def export(fig,name):
    for ext in ['svg','pdf','png']:fig.savefig(OUT/f'{name}.{ext}',dpi=301)
    fig.set_dpi(601);fig.canvas.draw()
    Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB').save(OUT/f'{name}.tiff',compression='tiff_lzw',dpi=(601,601))
    plt.close(fig)
def donorplot(ax,frame,categories,value,labels=None):
    for i,cat in enumerate(categories):
        for gi,g in enumerate(COL):
            d=frame[(frame.category==cat)&(frame.group==g)].sort_values('donor');y=d[value].to_numpy();off=(-.17,.17)[gi]
            ax.scatter(i+off+np.linspace(-.07,.07,len(y)),y,s=20,marker='o' if g=='Control' else '^',color=COL[g],edgecolor='white',linewidth=.35,zorder=3)
            if len(y):ax.plot([i+off-.10,i+off+.10],[y.mean()]*2,color=COL[g],lw=1.8,zorder=4)
    ax.set_xticks(range(len(categories)),[(labels or {}).get(x,x) for x in categories]);ax.set_xlim(-.55,len(categories)-.45);ax.grid(axis='y',alpha=.12)

# Figure 1 places sampling coverage and the confounded clinical comparison first.
fig=plt.figure(figsize=(174/25.4,186/25.4),layout='constrained')
gs=fig.add_gridspec(2,2,height_ratios=[1.35,1],width_ratios=[1.03,1],hspace=.09,wspace=.08)
ax=fig.add_subplot(gs[0,0]);cnt=pd.read_csv(R/'donor_celltype_qc.tsv',sep='\t').pivot(index='donor',columns='celltype',values='n_nuclei').reindex(index=DON,columns=CT)
ax.imshow(np.log10(cnt.clip(lower=1)),cmap='Blues',vmin=0,vmax=3.5,aspect='auto')
for i in range(12):
    for j in range(5):
        n=int(cnt.iloc[i,j]);ax.text(j,i,str(n),ha='center',va='center',fontsize=8,color='white' if n>500 else '#17222B')
        if n<20:ax.add_patch(Rectangle((j-.48,i-.48),.96,.96,fill=False,edgecolor='#A94F29',linewidth=1.1))
ax.set_xticks(range(5),CT);ax.set_yticks(range(12),[d.replace('Control','C').replace('UUO','U') for d in DON]);ax.axhline(6.5,color='white',lw=2)
panel(ax,'a','Donor coverage across lineages')
ax.text(0,-.14,'Outlined cells: fewer than 20 nuclei',transform=ax.transAxes,fontsize=8)
ax=fig.add_subplot(gs[0,1]);clinical=pd.read_csv(R/'donor_clinical_source_table.tsv',sep='\t')
for g in COL:
    d=clinical[clinical.group==g];ax.scatter(d.baseline_egfr_ml_min_1_73m2,d.tubulointerstitial_fibrosis_percent,marker='o' if g=='Control' else '^',color=COL[g],s=25,zorder=3)
    for _,z in d.iterrows():
        offset={'Control1':(3,-10),'Control6':(3,3),'Control2':(-16,-10),'Control4':(3,3),'Control5':(-14,3),'UUO2':(-16,4),'UUO5':(-16,-11),'UUO1':(-13,4)}.get(z.donor,(3,4))
        ax.annotate(z.donor.replace('Control','C').replace('UUO','U'),(z.baseline_egfr_ml_min_1_73m2,z.tubulointerstitial_fibrosis_percent),xytext=offset,textcoords='offset points',fontsize=8)
ax.set_xlabel('Baseline eGFR (mL/min/1.73 m²)');ax.set_ylabel('Tubulointerstitial fibrosis (%)');ax.set_ylim(-3,86);ax.set_xlim(33,106);ax.grid(alpha=.12)
panel(ax,'b','Clinical and histological context');leg(ax,loc='center right')
ax=fig.add_subplot(gs[1,0]);co=pd.read_csv(R/'captured_composition_comparisons.tsv',sep='\t')
for i,(ct,den) in enumerate([('TAL','fraction_all_captured'),('TAL','fraction_captured_epithelia'),('DCT','fraction_all_captured'),('DCT','fraction_captured_epithelia')]):
    z=co[(co.celltype==ct)&(co.denominator==den)].iloc[0];color='#446D87' if den=='fraction_all_captured' else '#777777'
    ax.plot([z.ci_low,z.ci_high],[i,i],color=color,lw=1.1);ax.scatter(z.effect,i,marker='o' if den=='fraction_all_captured' else 's',facecolor='white',edgecolor=color,s=25,zorder=3)
    ax.text(1.03,i,f'{z.fdr_permutation:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
ax.text(1.03,1.015,'FDR',transform=ax.transAxes,fontsize=8);ax.axvline(0,color='#AAA',ls=':',lw=.8)
ax.set_yticks(range(4),['TAL all','TAL epithelial','DCT all','DCT epithelial']);ax.set_ylim(3.7,-.8);ax.set_xlabel('Difference in capture fraction (pp)');panel(ax,'c','Composition depends on denominator')
ax=fig.add_subplot(gs[1,1]);st=pd.read_csv(R/'author_state_composition.tsv',sep='\t');st=st[st.celltype.isin(['TAL','DCT'])];ss=st[st.state.str.contains('Injured|Inflammatory')].groupby(['donor','group','celltype']).fraction_within_celltype.sum().rename('fraction').reset_index();ss=st[['donor','group','celltype']].drop_duplicates().merge(ss,how='left').fillna({'fraction':0});ss['percent']=ss.fraction*100;ss['category']=ss.celltype
donorplot(ax,ss,['TAL','DCT'],'percent');ax.set_ylim(-3,103);ax.set_ylabel('Injured / inflammatory nuclei (%)');panel(ax,'d','Author-labelled injury states')
export(fig,'Fig1');source('Fig1',['donor_celltype_qc.tsv','donor_clinical_source_table.tsv','captured_composition_comparisons.tsv','author_state_composition.tsv'],'The 12-donor comparison has unequal lineage coverage and a strong clinical/state contrast; capture estimates depend on their denominator.')

# Figure 2 gives all five primary estimates equal inferential visibility.
fig=plt.figure(figsize=(174/25.4,182/25.4),layout='constrained');gs=fig.add_gridspec(3,1,height_ratios=[1.1,1.2,1],hspace=.12)
top=gs[0].subgridspec(1,2,width_ratios=[1.3,1],wspace=.03)
ax=fig.add_subplot(top[0]);stats_ax=fig.add_subplot(top[1]);stats_ax.set_axis_off();stats_ax.set_xlim(0,1);stats_ax.set_ylim(4.6,-.7)
df=cmp.set_index('module').loc[PRIMARY]
for i,(m,z) in enumerate(df.iterrows()):
    if i%2==0:ax.axhspan(i-.45,i+.45,color='#F1F4F6',zorder=0)
    ax.plot([z.ci_low,z.ci_high],[i,i],color='#446D87',lw=1.3);ax.scatter(z.effect,i,s=28,facecolor='#446D87' if z.fdr<.05 else 'white',edgecolor='#446D87',zorder=3)
    for xpos,txt,ha in [(.015,f'{int(z.n_control)}/{int(z.n_uuo)}','left'),(.22,f'{z.effect:.2f} [{z.ci_low:.2f}, {z.ci_high:.2f}]','left'),(.985,f'{z.fdr:.3f}','right')]:
        stats_ax.text(xpos,i,txt,transform=stats_ax.get_yaxis_transform(),va='center',ha=ha,fontsize=8)
for xpos,txt,ha in [(.015,'n C/U','left'),(.22,'Difference [95% CI]','left'),(.985,'FDR','right')]:stats_ax.text(xpos,1.025,txt,transform=stats_ax.transAxes,ha=ha,fontsize=8)
ax.set_yticks(range(5),[SHORT[m] for m in PRIMARY]);ax.set_ylim(4.6,-.7);ax.set_xlim(-3,.99);ax.axvline(0,color='#AAA',ls=':',lw=.8);ax.set_xlabel('Mean log2 CPM difference');panel(ax,'a','Three primary programs have lower scores')
ax=fig.add_subplot(gs[1]);d=main[main.module.isin(PRIMARY)].copy();d['category']=d.module;donorplot(ax,d,PRIMARY,'score_mean_log2cpm',SHORT);ax.set_ylabel('Mean log2 CPM score');panel(ax,'b','Individual donor scores');leg(ax,'lower left',2)
ax=fig.add_subplot(gs[2]);lodo=pd.read_csv(R/'module_lodo_summary.tsv',sep='\t').set_index('module')
for i,m in enumerate(PRIMARY):
    z=lodo.loc[m];ax.plot([z.lodo_min,z.lodo_max],[i,i],color='#95AEBE',lw=3);ax.scatter(z.main_effect,i,color='#446D87',s=22,zorder=3)
    ax.text(.985,i,f'{int(z.n_same_direction)}/{int(z.n_lodo)}',transform=ax.get_yaxis_transform(),va='center',ha='right',fontsize=8)
ax.set_yticks(range(5),[SHORT[m] for m in PRIMARY]);ax.set_ylim(4.6,-.6);ax.set_xlim(-2.0,.45);ax.axvline(0,color='#AAA',ls=':',lw=.8);ax.set_xlabel('Mean log2 CPM difference');panel(ax,'c','Donor-omission ranges');ax.text(.985,1.04,'Same sign',transform=ax.transAxes,ha='right',fontsize=8)
export(fig,'Fig2');source('Fig2',['module_comparisons_primary_run.tsv','module_scores_all_scenarios.tsv','module_lodo_summary.tsv'],'All five prespecified programs are displayed; three meet program FDR, and donor omission describes stability rather than replication.')

# Figure 3 reports every fixed member in its own primary lineage.
fig,axes=plt.subplots(1,2,figsize=(174/25.4,142/25.4),layout='constrained',gridspec_kw={'wspace':.16})
genes=pd.read_csv(R/'fixed_transport_gene_results.tsv',sep='\t');mods=json.loads((ROOT/'03_protocol/modules_v1.json').read_text())
for j,(ax,ct) in enumerate(zip(axes,['TAL','DCT'])):
    members=list(dict.fromkeys(g for mm in mods if mm['celltype']==ct and mm['family']=='primary_transport' for g in mm['genes']))
    for i,g in enumerate(members):
        z=genes[(genes.celltype==ct)&(genes.gene==g)]
        if i%2==0:ax.axhspan(i-.48,i+.48,color='#F4F6F7',zorder=0)
        if len(z):
            z=z.iloc[0];ax.plot([z['CI.L'],z['CI.R']],[i,i],color='#536F81',lw=1.1);ax.scatter(z.logFC,i,s=24,facecolor='#446D87' if z.FDR_family<.05 else 'white',edgecolor='#446D87',zorder=3)
            if g in ['SLC12A1','SLC12A3','TRPM6']:ax.text(.98,i,f'q={z.FDR_family:.3f}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=8)
        else:ax.text(.05,i,'Below gene filter',transform=ax.get_yaxis_transform(),va='center',fontsize=8,color='#666')
    ax.set_yticks(range(len(members)),members,fontstyle='italic');ax.set_ylim(len(members)-.7,-.7);ax.axvline(0,color='#AAA',ls=':',lw=.8);ax.set_xlabel('Gene log2 fold change');ax.set_xlim(-7.2,3.5);panel(ax,'ab'[j],ct+' fixed transport genes')
export(fig,'Fig3');source('Fig3',['fixed_transport_gene_results.tsv','../03_protocol/modules_v1.json'],'TRPM6 meets the joint gene threshold; SLC12A1 and SLC12A3 have negative estimates with q about 0.057.')

fig,axes=plt.subplots(1,2,figsize=(174/25.4,106/25.4),layout='constrained',gridspec_kw={'wspace':.12})
order=['PC_water','PC_NaK','ICA_acid','ICB_base'];ax=axes[0]
for i,m in enumerate(order):
    z=cmp[cmp.module==m].iloc[0]
    if pd.notna(z.ci_low):
        ax.plot([z.ci_low,z.ci_high],[i,i],color='#536F81');ax.scatter(z.effect,i,facecolor='white',edgecolor='#446D87',s=25);txt=f'{z.fdr:.3f}'
    else:ax.scatter(z.effect,i,marker='x',color='#777',s=25);txt='NA'
    ax.text(1.03,i,txt,transform=ax.get_yaxis_transform(),va='center',fontsize=8)
ax.text(1.03,1.03,'FDR',transform=ax.transAxes,fontsize=8);ax.set_yticks(range(4),[SHORT[m] for m in order]);ax.set_ylim(3.7,-.7);ax.axvline(0,color='#AAA',ls=':',lw=.8);ax.set_xlabel('Mean log2 CPM difference');panel(ax,'a','Collecting-duct estimates')
ax=axes[1];d=main[main.module.isin(order)].copy();d['category']=d.module;donorplot(ax,d,order,'score_mean_log2cpm',SHORT);ax.tick_params(axis='x',labelrotation=35);ax.set_ylabel('Mean log2 CPM score');panel(ax,'b','Collecting-duct donor scores');leg(ax,'lower right')
export(fig,'FigS4');source('FigS4',['module_comparisons_primary_run.tsv','module_scores_all_scenarios.tsv'],'Collecting-duct estimates remain exploratory and uncertain, and IC-B is descriptive.')
(DEV/'00_admin/figure_contract.json').write_text(json.dumps(dict(width_mm=174,text_minimum_pt=8,raster_dpi=601,export_formats=['svg','pdf','png','tiff'],primary_results_changed=False,figures=data_manifest,preserved_supplementary_figures=['FigS1','FigS2','FigS3']),indent=2))
print('Exported Fig1 Fig2 Fig3 and FigS4')
