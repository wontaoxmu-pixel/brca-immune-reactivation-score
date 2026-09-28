from pathlib import Path
import csv,sys,json,numpy as np,pandas as pd
root=Path(__file__).resolve().parents[2]; src=root/'data/processed/scRNA';out=root/'results/revision_R1/tables'
sets={'CD8':['CD8A','CD8B','PDCD1','LAG3','HAVCR2','TIGIT','TOX','CXCL13','GZMB','PRF1','NKG7'],'AP':['HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','TAPBP','PSMB8','PSMB9','NLRC5'],'IFN':['IFNG','STAT1','IRF1','CXCL9','CXCL10','GBP1','GBP5','ISG15','IFIT1','MX1']}
markers=['CD3D','CD3E','TRAC','CD4','CD8A','CD8B','NKG7','GNLY','FOXP3','IL7R','CCR7','TCF7','GZMK','GZMB','ZFP36','LAG3','PDCD1','HAVCR2','TOX','IFNG','IFIT1','CD68','LYZ','LST1','FCGR3A','S100A9','CD1C','CLEC9A','EPCAM','KRT8','KRT18','PTPRC']
genes=[s.strip().upper() for s in (src/'count_matrix_genes.tsv').open()];barcodes=[s.strip() for s in (src/'count_matrix_barcodes.tsv').open()]
m=pd.read_csv(src/'metadata.csv');assert list(m.iloc[:,0])==barcodes
m['sample']=m['orig.ident'];m['group']=m['sample']+'|'+m.celltype_subset
sel=m.subtype.eq('TNBC');groups=sorted(m.loc[sel,'group'].unique());idx={g:i for i,g in enumerate(groups)}
cellgroup=np.array([idx.get(g,-1) if keep else -1 for g,keep in zip(m.group,sel)])
wanted=list(dict.fromkeys(sum(sets.values(),[])+markers));present=[g for g in wanted if g in genes];gmap={str(genes.index(g)+1).encode():i for i,g in enumerate(present)}
cache=root/'data/revision_R1/scRNA_selected_counts.npz'
if cache.exists():
 a=np.load(cache);counts=a['counts'];det=a['det'];total=a['total'];nnz=int(a['nnz'])
else:
 counts=np.zeros((len(groups),len(present)));det=np.zeros_like(counts);total=np.zeros(len(m));nnz=0
 with (src/'count_matrix_sparse.mtx').open('rb') as f:
  for line in f:
   if not line.startswith(b'%'):nr,nc,declared=map(int,line.split());break
  assert nr==len(genes) and nc==len(m)
  for line in f:
   a,b,v=line.split();j=int(b)-1;v=float(v);total[j]+=v;nnz+=1;k=cellgroup[j]
   if k>=0 and a in gmap:counts[k,gmap[a]]+=v;det[k,gmap[a]]+=1
   if nnz%30000000==0:print(nnz,flush=True)
 assert nnz==declared
 np.savez_compressed(cache,counts=counts,det=det,total=total,nnz=nnz)
# Verify supplied UMI totals rather than assume same assay.
qc={'cells':len(m),'nonzero_entries':nnz,'max_UMI_difference':float(np.max(np.abs(total-m.nCount_RNA.to_numpy()))),'TNBC_identifiers':int(m.loc[sel,'sample'].nunique()),'patient_identity_status':'GEO SOFT sample titles verified against source tumor identifiers; sample-level biological units'}
(out/'scRNA_input_qc.json').write_text(json.dumps(qc,indent=2))
m['actual_total']=total;meta=m.loc[sel].groupby('group').agg(sample=('sample','first'),subset=('celltype_subset','first'),major=('celltype_major','first'),n_cells=('sample','size'),total_umi=('actual_total','sum')).loc[groups].reset_index()
cpm=counts/meta.total_umi.to_numpy()[:,None]*1e6; log=np.log2(cpm+1)
for k,g in sets.items():meta[k]=log[:,[present.index(x) for x in g]].mean(1)
meta['composite']=meta[list(sets)].mean(1);meta.to_csv(out/'scRNA_sample_subset_scores.csv',index=False)
rows=[]
for i,r in meta.iterrows():
 for j,g in enumerate(present):rows.append(dict(sample=r['sample'],subset=r['subset'],major=r['major'],n_cells=r.n_cells,gene=g,counts=counts[i,j],cpm=cpm[i,j],log2cpm=log[i,j],detection=det[i,j]/r.n_cells))
pd.DataFrame(rows).to_csv(out/'scRNA_marker_expression.csv',index=False)
rows=[]
for th in [10,20,50]:
 for state,d in meta[meta.n_cells>=th].groupby('subset'):
  for k in list(sets)+['composite']:rows.append(dict(threshold=th,subset=state,component=k,n_samples=d['sample'].nunique(),median=d[k].median(),min=d[k].min(),max=d[k].max()))
pd.DataFrame(rows).to_csv(out/'scRNA_threshold_sensitivity.csv',index=False)
# Simulated RNA-contribution mixtures: fixed source profiles, no resampling of cells as biological replicates.
sim=[]
for sample,indices in meta.groupby('sample').groups.items():
 ids=list(indices);ok=[i for i in ids if meta.loc[i,'n_cells']>=20]
 def aggregate(ids):
  return counts[ids].sum(0)/meta.loc[ids,'total_umi'].sum()*1e6 if ids else None
 mal=[i for i in ok if meta.loc[i,'major']=='Cancer Epithelial']; ts=[i for i in ok if 'CD8+' in meta.loc[i,'subset']]; my=[i for i in ok if meta.loc[i,'major']=='Myeloid']
 profiles=[aggregate(x) for x in [mal,ts,my]]
 if all(x is not None for x in profiles):
  a,t,y=profiles
  for f in np.linspace(0,0.8,17):
   for scenario,mix in [('CD8_RNA_fraction',(1-f)*a+f*t),('myeloid_RNA_fraction',(1-f)*a+f*y)]:
    vals={k:float(np.log2(mix[[present.index(g) for g in gs]]+1).mean()) for k,gs in sets.items()};sim.append(dict(sample=sample,scenario=scenario,fraction=f,**vals,composite=np.mean(list(vals.values()))))
  low=[i for i in ts if 'ZFP36' in meta.loc[i,'subset']];high=[i for i in ts if 'LAG3' in meta.loc[i,'subset']]
  if low and high:
   lo,hi=aggregate(low),aggregate(high)
   for f in np.linspace(0,1,21):
    mix=.7*a+.1*y+.2*((1-f)*lo+f*hi);vals={k:float(np.log2(mix[[present.index(g) for g in gs]]+1).mean()) for k,gs in sets.items()};sim.append(dict(sample=sample,scenario='LAG3_share_fixed_CD8_RNA_20pct',fraction=f,**vals,composite=np.mean(list(vals.values()))))
simdf=pd.DataFrame(sim)
# Fixed standardization references are the eligible sample-by-subset profiles, not varied across mixture points.
mu=log[meta.n_cells.to_numpy()>=20].mean(0); sd=log[meta.n_cells.to_numpy()>=20].std(0,ddof=1)
# Component log-CPM simulations are primary descriptive outputs; do not label them the cohort-z standardized bulk score.
simdf['score_scale']='mean component log2(CPM+1); not bulk cohort z-score'
simdf.to_csv(out/'scRNA_composition_simulations.csv',index=False)
print(qc);print('mixture rows',len(sim));print('finished',flush=True)
