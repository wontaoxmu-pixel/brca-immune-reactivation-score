from pathlib import Path
import json,hashlib,gzip,re
import pandas as pd,numpy as np
root=Path(__file__).resolve().parents[2];t=root/'results/revision_R1/tables';audit=root/'manuscript/plos_one_revision_R1/audit'
flows=[];cors=[];vifs=[]
for ds in ['GSE58812','GSE96058','METABRIC','TCGA']:
 a=pd.read_csv(t/f'{ds}_analysis_rows.csv');assert not a.id.duplicated().any();assert set(a.event.dropna()).issubset({0,1})
 ok=a[['time','event','age','original','immune']].notna().all(axis=1)&a.time.gt(0)&a.event.isin([0,1])
 flows.append(dict(dataset=ds,expression_samples=len(a),duplicate_ids=int(a.id.duplicated().sum()),missing_age=int(a.age.isna().sum()),missing_time=int(a.time.isna().sum()),missing_event=int(a.event.isna().sum()),nonpositive_time=int(a.time.le(0).sum()),eligible=int(ok.sum()),events=int(a.loc[ok,'event'].sum()),stage_observed=int(a.stage.notna().sum()),purity_observed=int(a.purity.notna().sum()),median_observed_time_days=float(a.loc[ok,'time'].median())))
 for sc in ['gene_equal','singscore']:
  cors.append(dict(dataset=ds,score_a='original',score_b=sc,spearman=a.original.corr(a[sc],method='spearman')))
 for ss,b in [('all',a[ok]),('tnbc',a[ok & a.tnbc.eq(True)]),('basal',a[ok & a.basal.eq(True)])]:
  if len(b)<20:continue
  x=np.column_stack([np.ones(len(b)),b.age,b.immune]);y=b.original.to_numpy();fit=x@np.linalg.lstsq(x,y,rcond=None)[0];r2=1-np.sum((y-fit)**2)/np.sum((y-y.mean())**2)
  vifs.append(dict(dataset=ds,subset=ss,n=len(b),score_VIF_with_age_and_immune=1/(1-r2)))
pd.DataFrame(flows).to_csv(t/'sample_flow_and_covariates.csv',index=False);pd.DataFrame(cors).to_csv(t/'scoring_method_correlations.csv',index=False);pd.DataFrame(vifs).to_csv(t/'collinearity_audit.csv',index=False)
soft=gzip.open(root/'data/revision_R1/GSE176078_family.soft.gz','rt').read();pairs=[]
for block in soft.split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0].strip();m=re.search(r'!Sample_title = (.+)',block)
 if m:pairs.append({'GSM':gsm,'source_title':m.group(1)})
pd.DataFrame(pairs).to_csv(t/'GSE176078_source_sample_titles.csv',index=False)
qc=json.loads((t/'scRNA_input_qc.json').read_text());qc['sample_mapping_evidence']='GSE176078_family.soft.gz sample titles; source article describes 26 primary tumors including 10 TNBC. Analyses use sample identifiers; no assertion of independent repeated specimens.';(t/'scRNA_input_qc.json').write_text(json.dumps(qc,indent=2))
print(pd.DataFrame(flows).to_string(index=False));print(pd.DataFrame(cors).to_string(index=False))
