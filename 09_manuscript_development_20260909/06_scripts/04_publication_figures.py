#!/usr/bin/env python3
"""Reuse the validated numerical plots with publication sizing and labels."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.lines import Line2D
from PIL import Image
DEV=Path(__file__).resolve().parents[1];ROOT=DEV.parent;OUT=DEV/'05_figures'
NAME={'Figure1_donor_QC':'FigS1','Figure2_primary_transport':'Fig1','Figure3_genes_collecting_duct':'Fig2','Figure4_stability_covariation':'FigS2'}
TITLES={'Donor contribution by library (%)':'Library contribution (%)','Nuclei per donor and lineage':'Nuclei by donor and lineage','Classical markers in reference donors':'Reference lineage markers','RNA quality across donors':'RNA counts by donor','Relative capture composition':'Capture composition','Author-labelled injury states':'Author-labelled injury states','Prespecified transport effects':'Primary transport programs','Individual donor expression':'Donor program scores','TAL: prespecified member genes':'TAL member genes','DCT: prespecified member genes':'DCT member genes','Collecting-duct program estimates':'Collecting duct estimates','Collecting-duct donor expression':'Collecting duct scores','Leave-one-donor-out estimates':'Donor omission','Prespecified sensitivity estimates':'Sensitivity estimates','TAL transport and epithelial injury':'TAL salt and injury','DCT transport and epithelial injury':'DCT NaCl and injury'}
REPLACE={'Gene log2 fold change (UUO - reference)':'Gene log2 fold change','Difference in mean log2 CPM':'Mean log2 CPM difference','Epithelial injury score (mean log2 CPM)':'Injury score (mean log2 CPM)','TAL salt score (mean log2 CPM)':'TAL salt score','DCT NaCl score (mean log2 CPM)':'DCT NaCl score','Filled: FDR <0.05 across TAL + DCT gene tests':'Filled: joint gene FDR <0.05','Injured + inflammatory states; descriptive annotation':'Injured/inflammatory; descriptive','TAL: 7 reference / 5 UUO; DCT: 6 reference / 5 UUO':'Reference / UUO: TAL 7/5; DCT 6/5','Bold counts: below the 20-nucleus threshold':'Bold: fewer than 20 nuclei','Size: donor-mean detection fraction (0-100%)':'Size: detection fraction (0-100%)','Bar: estimate range; dot: full-sample estimate':'Bar: range; dot: main estimate','Healthy-state DCT: only 1 eligible UUO donor':'Healthy DCT: 1 eligible UUO donor','Median RNA UMIs per nucleus':'Median RNA counts per nucleus','Module mean log2 CPM':'Mean log2 CPM score'}
def export(fig,name):
    short=NAME.get(name,name)
    height={'Fig1':184,'Fig2':200,'FigS1':196,'FigS2':200,'FigS3':116}[short]
    fig.set_size_inches(174/25.4,height/25.4)
    for item in fig.findobj(Text):
        item.set_fontsize(max(8,item.get_fontsize()))
        item.set_text(REPLACE.get(item.get_text(),item.get_text()))
    for ax in fig.axes:
        if ax.get_title(loc='left') in TITLES:ax.set_title(TITLES[ax.get_title(loc='left')],loc='left',fontsize=9,pad=10)
        if ax.get_legend() is not None:
            ax.get_legend().remove()
            ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color='#446D87',label='RCC reference',markersize=4),Line2D([],[],marker='^',linestyle='',color='#A94F29',label='UTUC + obstruction',markersize=4)],loc='best',fontsize=8)
    fig.savefig(OUT/f'{short}.svg',dpi=301)
    fig.savefig(OUT/f'{short}.pdf',dpi=301)
    fig.savefig(OUT/f'{short}.png',dpi=301)
    fig.set_dpi(601);fig.canvas.draw()
    im=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')
    im.save(OUT/f'{short}.tiff',compression='tiff_lzw',dpi=(601,601))
    plt.close(fig)

source=(ROOT/'04_scripts/06_figures.py').read_text()
source=source.replace("O=ROOT/'06_figures'", "O=DEV/'05_figures'")
source=source.replace("R/'marker_audit_donor_balanced.tsv'", "DEV/'02_analysis/marker_donor_balanced.tsv'")
start=source.index('def save(fig,name):');end=source.index('\ndef donorplot',start)
source=source[:start]+"def save(fig,name):\n    publication_export(fig,name)\n"+source[end:]
source=source.replace("'#547D99'","'#446D87'").replace("'#C76A43'","'#A94F29'")
source=source.replace("ax.scatter(i+off+j,y,s=15,color=COL[g]", "ax.scatter(i+off+j,y,s=15,marker='o' if g=='Control' else '^',color=COL[g]")
source=source.replace("ax.scatter(s.epithelial_injury,s[mo],color=COL[g]", "ax.scatter(s.epithelial_injury,s[mo],marker='o' if g=='Control' else '^',color=COL[g]")
source=source.replace("ax.text(0,-.46,", "ax.text(0,-.35,")
exec(compile(source,str(ROOT/'04_scripts/06_figures.py'),'exec'),{'__file__':str(ROOT/'04_scripts/06_figures.py'),'DEV':DEV,'publication_export':export})

data=pd.read_csv(DEV/'02_analysis/member_omission_estimates.tsv',sep='\t')
summary=pd.read_csv(DEV/'02_analysis/member_influence_summary.tsv',sep='\t').set_index('module')
fig,axes=plt.subplots(1,3,figsize=(174/25.4,116/25.4),layout='constrained')
for k,(ax,mod,label) in enumerate(zip(axes,['TAL_salt','DCT_NaCl','DCT_Mg'],['TAL salt','DCT NaCl','DCT magnesium'])):
    z=data[data.module.eq(mod)].reset_index(drop=True)
    for i,r in z.iterrows():
        ax.plot([r.ci_low,r.ci_high],[i,i],color='#677D8C',lw=1)
        ax.scatter(r.effect,i,color='#446D87',s=16,zorder=3)
    ax.axvline(0,color='#AAA',ls=':',lw=.7)
    ax.axvline(summary.loc[mod,'main_effect'],color='#777',ls='--',lw=.7)
    ax.set_yticks(range(len(z)),z.omitted_gene);ax.set_ylim(len(z)-.5,-.5)
    ax.set_xlabel('Mean log2 CPM difference');ax.set_title(label,loc='left',fontsize=9)
    ax.text(-.22,1.05,'abc'[k],transform=ax.transAxes,fontweight='bold',fontsize=11)
axes[0].set_ylabel('Omitted member')
export(fig,'FigS3')
(DEV/'07_qa/figure_generation_receipt.json').write_text(json.dumps({'figures':list(NAME.values())+['FigS3'],'width_mm':174,'raster_export_dpi':601,'backend':'matplotlib','authoritative_data':'unchanged first-round result tables; new post hoc influence estimates','primary_statistics_changed':False},indent=2))
print('Generated five publication figures')
