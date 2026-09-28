from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullFormatter
from matplotlib.lines import Line2D
from PIL import Image
root=Path(__file__).resolve().parents[2];t=root/'results/revision_R1/tables';out=root/'results/revision_R1/figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
colors={'age':'#737373','score':'#2878B5','immune':'#C67B21','both':'#8F3985'}
labels={'tnbc':'TNBC','pathology_tnbc':'Pathology TNBC','pam50_basal':'Basal-like','all_primary':'All primary','tnbc_MFS':'TNBC (MFS)'}
def lab(ds,ss):return ds+' / '+labels.get(ss,ss)
def save(fig,name):
 fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(out/(name+'.tif'),dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'});plt.close(fig); im=Image.open(out/(name+'.tif')).convert('RGB'); im.save(out/(name+'.tif'),compression='tiff_lzw',dpi=(600,600))
def panel(ax,s):ax.text(-.06,1.06,s,transform=ax.transAxes,weight='bold',fontsize=11,va='top')
co=pd.read_csv(t/'abundance_correlations.csv');cx=pd.read_csv(t/'cox_results.csv');bt=pd.read_csv(t/'bootstrap_performance.csv');inc=pd.read_csv(t/'incremental_results.csv')
fig=plt.figure(figsize=(7.1,5.3));gs=GridSpec(2,1,height_ratios=[1,2],hspace=.6);ax=fig.add_subplot(gs[0]);ax.axis('off');panel(ax,'A')
boxes=[(.13,'Fixed 31 genes\n3 biological components'),(.5,'Bulk cohorts\nAbundance + survival'),(.87,'Single-cell samples\nStates + RNA mixtures')]
for x,s in boxes:ax.text(x,.57,s,ha='center',va='center',bbox=dict(boxstyle='round,pad=.65',fc='#EEF3F8',ec='#8FA6BA'),fontsize=9)
for a,b in [(.28,.34),(.66,.72)]:ax.annotate('',xy=(b,.57),xytext=(a,.57),arrowprops=dict(arrowstyle='->',color='#586B7E'))
ax.text(.5,-.04,'Separate association, incremental information and biological interpretation',ha='center',fontsize=9)
ax=fig.add_subplot(gs[1]);panel(ax,'B');ab=['immune','T cells','CD8 T cells','Cytotoxic lymphocytes','Monocytic lineage'];ds=['GSE58812','GSE96058','METABRIC','TCGA'];a=co[(co.score=='original') & co.abundance.isin(ab)].pivot(index='dataset',columns='abundance',values='rho').reindex(index=ds,columns=ab)
im=ax.imshow(a,vmin=0,vmax=1,cmap='YlGnBu',aspect='auto');ax.set_yticks(range(4),ds);ax.set_xticks(range(5),['ESTIMATE\nImmuneScore','MCP\nT cells','MCP\nCD8 T cells','MCP\nCytotoxic','MCP\nMonocytic']);ax.set_title('Correlation with expression-based abundance summaries',pad=10)
for i in range(4):
 for j in range(5):ax.text(j,i,f'{a.iloc[i,j]:.2f}',ha='center',va='center',color='white' if a.iloc[i,j]>.65 else 'black')
fig.colorbar(im,ax=ax,pad=.02,label="Spearman's rho");save(fig,'Fig1')
# Ordered identical forest rows, one parent-cohort SD.
a=cx[(cx.score=='original')&(cx.adjustment=='age')].copy();a=a[~a.subset.eq('tnbc_MFS')];order=[(d,s) for d in ds for s in ['all_primary','tnbc','pathology_tnbc','pam50_basal'] if ((a.dataset==d)&(a.subset==s)).any()]
fig,axes=plt.subplots(1,2,figsize=(7.1,5.7),sharey=True,gridspec_kw={'wspace':.13})
for ax,adj,title,pp in zip(axes,['age','immune'],['Age adjusted','Age + ImmuneScore adjusted'],['A','B']):
 panel(ax,pp);v=cx[(cx.score=='original')&(cx.adjustment==adj)].set_index(['dataset','subset']);
 for i,k in enumerate(order):
  r=v.loc[k];color='#2878B5' if k[1]!='all_primary' else '#7F8C8D';ax.errorbar(r.hr,i,xerr=[[r.hr-r.lower],[r.upper-r.hr]],fmt='o',ms=4,color=color,capsize=2)
 ax.axvline(1,color='#888',ls='--',lw=.8);ax.set_xscale('log');ax.xaxis.set_minor_formatter(NullFormatter());ax.set_xlim(.18,3.3);ax.set_xticks([.25,.5,1,2,3],[.25,.5,1,2,3]);ax.set_title(title,fontsize=9);ax.set_xlabel('HR per parent-cohort SD (95% CI)');ax.grid(axis='x',alpha=.15)
axes[0].set_yticks(range(len(order)),[lab(*k) for k in order]);axes[0].invert_yaxis();fig.tight_layout();save(fig,'Fig2')
# Performance all analysis sets, panels retain matched rows.
order=[(d,s) for d in ds for s in ['all_primary','tnbc','pathology_tnbc','pam50_basal','tnbc_MFS'] if ((bt.dataset==d)&(bt.subset==s)).any()]
fig,axes=plt.subplots(1,2,figsize=(7.1,6.2),sharey=True,gridspec_kw={'width_ratios':[1.3,1],'wspace':.18})
for idx,k in enumerate(order):
 for off,m in zip([-.20,0,.20],['score','immune','both']):
  r=bt[(bt.dataset==k[0])&(bt.subset==k[1])&(bt.metric==m)].iloc[0];axes[0].hlines(idx+off,r.lower,r.upper,color=colors[m],lw=.8);axes[0].plot(r.estimate,idx+off,'o',ms=3,color=colors[m],label=m if idx==0 else None)
 r=bt[(bt.dataset==k[0])&(bt.subset==k[1])&(bt.metric=='delta_both_minus_immune')].iloc[0];axes[1].hlines(idx,r.lower,r.upper,color=colors['both'],lw=1);axes[1].plot(r.estimate,idx,'o',color=colors['both'],ms=4);axes[1].plot(r.corrected,idx,'x',color='black',ms=4)
axes[0].set_yticks(range(len(order)),[lab(*k) for k in order]);axes[0].invert_yaxis();axes[0].set_xlabel("Apparent Harrell's C-index (95% CI)");axes[0].set_xlim(.4,.95);axes[0].legend(handles=[Line2D([0],[0],color=colors[m],marker='o',ms=3,label=l) for m,l in [('score','Score + age'),('immune','ImmuneScore + age'),('both','Both + age')]],loc='upper center',bbox_to_anchor=(.5,-.10),fontsize=8,ncol=1,frameon=False)
axes[1].axvline(0,color='#999',ls='--',lw=.8);axes[1].set_xlabel('Delta C-index vs ImmuneScore + age');axes[1].set_title('Circle: apparent; cross: optimism corrected',fontsize=8)
for ax,pp in zip(axes,['A','B']):panel(ax,pp);ax.grid(axis='x',alpha=.15)
fig.tight_layout();save(fig,'Fig3')
sc=pd.read_csv(t/'scRNA_sample_subset_scores.csv');sim=pd.read_csv(t/'scRNA_composition_simulations.csv');eligible=sc[sc.n_cells>=20]
states=['T_cells_c4_CD8+_ZFP36','T_cells_c5_CD8+_GZMK','T_cells_c7_CD8+_IFNG','T_cells_c8_CD8+_LAG3','Myeloid_c1_LAM1_FABP5','Myeloid_c9_Macrophage_2_CXCL10']
# Add original malignant subsets as an explicitly sample-aggregated group for context.
h=eligible[eligible.subset.isin(states)].groupby('subset')[['CD8','AP','IFN']].median().reindex(states);ns=eligible.groupby('subset')['sample'].nunique()
fig=plt.figure(figsize=(7.1,6.8));gs=GridSpec(2,2,hspace=.65,wspace=.4);ax=fig.add_subplot(gs[0,:]);panel(ax,'A');im=ax.imshow(h,aspect='auto',cmap='viridis');ax.set_xticks([0,1,2],['CD8-associated','Antigen presentation','IFN response']);ax.set_yticks(range(len(states)),[s.replace('T_cells_','').replace('Myeloid_','')+f' (n={ns.get(s,0)})' for s in states]);fig.colorbar(im,ax=ax,label='Median component log2(CPM + 1)',pad=.02)
for j,(scenario,title) in enumerate([('CD8_RNA_fraction','Changing lineage RNA contribution'),('LAG3_share_fixed_CD8_RNA_20pct','Changing CD8-state composition')]):
 ax=fig.add_subplot(gs[1,j]);panel(ax,['B','C'][j]);use=[scenario,'myeloid_RNA_fraction'] if j==0 else [scenario]
 for scenario,color in zip(use,['#2878B5','#C67B21']):
  a=sim[sim.scenario==scenario]
  for s,d in a.groupby('sample'):ax.plot(d.fraction,d.composite,color=color,alpha=.17,lw=.7)
  b=a.groupby('fraction').composite.median();ax.plot(b.index,b.values,color=color,lw=2,label=('CD8' if scenario=='CD8_RNA_fraction' else 'Myeloid' if scenario=='myeloid_RNA_fraction' else 'LAG3 vs ZFP36'))
 ax.set_xlabel('RNA contribution fraction' if j==0 else 'LAG3 share within fixed 20% CD8 RNA');ax.set_ylabel('Descriptive composite log-CPM score');ax.set_title(title,fontsize=8);ax.legend(fontsize=8)
fig.subplots_adjust(left=.31,right=.97,top=.95,bottom=.09);save(fig,'Fig4')
# Supplementary alternative scoring associations
fig,ax=plt.subplots(figsize=(7.1,4.8));panel(ax,'A');a=cx[(cx.adjustment=='age')&(~cx.subset.eq('tnbc_MFS'))];order=[(d,s) for d in ds for s in ['all_primary','tnbc','pathology_tnbc','pam50_basal'] if ((a.dataset==d)&(a.subset==s)).any()]
for k,(scname,color,off) in enumerate([('original','#2878B5',-.2),('gene_equal','#C67B21',0),('singscore','#8F3985',.2)]):
 for i,key in enumerate(order):
  r=a[(a.dataset==key[0])&(a.subset==key[1])&(a.score==scname)].iloc[0];ax.errorbar(r.hr,i+off,xerr=[[r.hr-r.lower],[r.upper-r.hr]],fmt='o',ms=3,color=color,label=scname if i==0 else None)
ax.axvline(1,color='#999',ls='--');ax.set_xscale('log');ax.xaxis.set_minor_formatter(NullFormatter());ax.set_xticks([.5,1,1.5],[.5,1,1.5]);ax.set_yticks(range(len(order)),[lab(*k) for k in order]);ax.invert_yaxis();ax.set_xlabel('Age-adjusted HR per parent-cohort SD (95% CI)');ax.legend(fontsize=8);fig.tight_layout();save(fig,'S1_Fig')
# S2 explicit sample x-axis; existing component scores.
a=pd.read_csv(root/'results/tables/GSE58812_immune_reactivation_scores.csv').sort_values('immune_reactivation_score');cols=['exhausted_cd8_t_cell','antigen_presentation','interferon_response'];fig,ax=plt.subplots(figsize=(7.1,2.4));im=ax.imshow(a[cols].to_numpy().T,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5);ax.set_yticks(range(3),['CD8-associated','Antigen presentation','IFN response']);ax.set_xlabel('GSE58812 samples ordered by original composite score');ax.set_xticks([0,26,53,80,106],[1,27,54,81,107]);fig.colorbar(im,ax=ax,label='Component z-score');fig.tight_layout();save(fig,'S2_Fig')
# Threshold counts by source subset, not cell-level significance.
th=pd.read_csv(t/'scRNA_threshold_sensitivity.csv');a=th[(th.component=='composite')&th.subset.isin(states)];fig,ax=plt.subplots(figsize=(7.1,3.4))
for i,s in enumerate(states):
 b=a[a.subset==s].set_index('threshold').reindex([10,20,50]);ax.plot([10,20,50],b.n_samples,marker='o',label=s.replace('T_cells_','').replace('Myeloid_',''))
ax.set_xlabel('Minimum cells per sample and subset');ax.set_ylabel('Eligible tumor samples');ax.set_xticks([10,20,50]);ax.set_ylim(0,11);ax.legend(fontsize=8,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.02),frameon=False);fig.tight_layout();save(fig,'S3_Fig')
params=[]
for p in out.glob('*.tif'):
 im=Image.open(p);params.append({'file':p.name,'width_px':im.width,'height_px':im.height,'dpi':im.info.get('dpi'),'mode':im.mode,'bytes':p.stat().st_size})
pd.DataFrame(params).to_csv(t/'figure_parameters.csv',index=False)
print('figures complete')
