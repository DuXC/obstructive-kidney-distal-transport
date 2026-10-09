#!/usr/bin/env python3
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parents[1]
RES=Path(sys.argv[1]) if len(sys.argv)>1 else P/'05_results'
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'], 'font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.7,'legend.frameon':False})
colors={1:'#356F85',5:'#C16B3A'};markers={1:'o',5:'s'}
gene=pd.read_csv(RES/'gene_primary_contrasts.tsv',sep='\t')
prog=pd.read_csv(RES/'program_primary_contrasts.tsv',sep='\t')
def export(fig,name):
    for ext in ['svg','pdf','png','tiff']:
        kw={'dpi':600 if ext=='tiff' else 300} if ext in ['png','tiff'] else {}
        if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(P/'06_figures'/(name+'.'+ext),**kw)
    plt.close(fig)

fig,axes=plt.subplots(2,3,figsize=(174/25.4,149/25.4))
fig.subplots_adjust(left=.11,right=.98,bottom=.15,top=.87,wspace=.52,hspace=.62)
items=[(gene,'SLC12A1','Slc12a1'),(gene,'SLC12A3','Slc12a3'),(gene,'TRPM6','Trpm6'),(prog,'TAL_CaMg','TAL Ca/Mg program'),(prog,'DCT_Mg','DCT Mg program'),(prog,'DCT_Ca','DCT Ca program')]
for i,(ax,(df,entity,label)) in enumerate(zip(axes.flat,items)):
    d=df[df.entity.eq(entity)&df.kind.eq('same_time')]
    ax.axhline(0,color='#777777',lw=.6,zorder=0)
    for duration in [1,5]:
        r=d[d.duration_days.eq(duration)].sort_values('recovery_days')
        xs=np.array([0,1,2])+(-.11 if duration==1 else .11)
        for x,(_,v) in zip(xs,r.iterrows()):
            ax.errorbar(x,v.effect,yerr=[[v.effect-v.ci_low],[v.ci_high-v.effect]],color=colors[duration],marker=markers[duration],mfc=colors[duration] if v.fdr<.05 else 'white',mec=colors[duration],markersize=3.8,lw=.9,elinewidth=.9,capsize=2,ls='none')
    ax.set_xticks([0,1,2],['Acute','10','28']);ax.set_xlim(-.45,2.45)
    ax.set_title(label,fontsize=8.5,pad=8,fontstyle='italic' if i<3 else 'normal')
    ax.text(-.23,1.11,chr(97+i),transform=ax.transAxes,fontweight='bold',fontsize=9)
    if i>=3:ax.set_xlabel('Recovery day',fontsize=8,labelpad=5)
    if i%3==0:ax.set_ylabel('Affected minus sham\n(normalized log2 signal)',fontsize=8,labelpad=5)
fig.legend([Line2D([0],[0],marker=markers[d],color=colors[d],ls='none',markersize=4) for d in [1,5]],['1-day obstruction history','5-day obstruction history'],loc='upper center',bbox_to_anchor=(.55,.995),ncol=2,fontsize=8)
fig.text(.11,.047,'Filled markers: FDR < 0.05; open markers: FDR ≥ 0.05.',fontsize=8)
fig.text(.11,.022,'Acute and recovery samples are different animals; colours identify obstruction duration.',fontsize=7.5)
export(fig,'Fig5')

fig=plt.figure(figsize=(174/25.4,158/25.4))
gs=fig.add_gridspec(2,2,left=.10,right=.98,top=.90,bottom=.10,wspace=.40,hspace=.52,height_ratios=[1,1.10])
labels=['1 d / 0 d','1 d / 10 d','1 d / 28 d','5 d / 0 d','5 d / 10 d','5 d / 28 d']
contrasts=[f'd{d}_r{r}_vs_sham' for d in [1,5] for r in [0,10,28]]
for j,mode in enumerate(['detectable','date_adjusted']):
    ax=fig.add_subplot(gs[0,j]);df=pd.read_csv(RES/f'program_{mode}_contrasts.tsv' if mode=='detectable' else RES/'program_date_adjusted.tsv',sep='\t')
    df=df[df.entity.eq('DCT_Mg')].set_index('contrast').reindex(contrasts)
    ax.axvline(0,color='#888888',lw=.65)
    for y,(_,r) in enumerate(df.iterrows()):
        color=colors[1 if y<3 else 5]
        if r.status=='tested':ax.errorbar(r.effect,y,xerr=[[r.effect-r.ci_low],[r.ci_high-r.effect]],fmt='o',markersize=3.5,color=color,lw=.8,capsize=2)
        else:ax.text(-.92,y,'Not evaluable',va='center',fontsize=7.5,color='#6C6C6C')
    ax.set_yticks(range(6),labels);ax.set_ylim(5.55,-.55);ax.set_xlim(-1.05,.30)
    ax.set_xlabel('DCT Mg score difference (95% CI)',fontsize=8)
    ax.set_title('Above-background sensitivity' if j==0 else 'Acquisition-date sensitivity',fontsize=8.5,pad=10)
    ax.text(-.26,1.10,chr(97+j),transform=ax.transAxes,fontweight='bold',fontsize=9)

ax=fig.add_subplot(gs[1,:]);g=pd.read_csv(RES/'gene_values.tsv',sep='\t')
d=g[g.analysis.eq('detectable')].groupby('mouse_symbol').value.count().sort_values()
ax.bar(np.arange(len(d)),d.values,color='#A3BBC5',width=.68)
ax.set_ylim(0,79);ax.set_yticks([0,20,40,60,71]);ax.axhline(71,color='#999999',lw=.6,ls=':')
ax.set_xticks(np.arange(len(d)),d.index,rotation=45,ha='right',fontsize=8,fontstyle='italic')
ax.set_ylabel('Arrays with detectable gene\n(out of 71)',fontsize=8)
for x,n in enumerate(d.values):ax.text(x,n+1,str(n),ha='center',va='bottom',fontsize=7.5)
ax.set_title('At least one clean probe spot above background',fontsize=8.5,pad=10)
ax.text(-.095,1.12,'c',transform=ax.transAxes,fontweight='bold',fontsize=9)
export(fig,'FigS7')
print('Exported Fig5 and FigS7 as SVG, PDF, TIFF and PNG')
