set.seed(123);.libPaths(c('.Rlib',.libPaths()));suppressPackageStartupMessages({library(data.table);library(survival)})
t='results/revision_R1/tables';expected=fread(file.path(t,'cox_results.csv'))[score=='original' & adjustment%in%c('age','immune')];checks=list()
for(ds in unique(expected$dataset)){
 d=fread(file.path(t,paste0(ds,'_analysis_rows.csv')))
 for(ss in unique(expected[dataset==ds]$subset)){
  a=copy(d);if(ss=='pathology_tnbc')a=a[tnbc==TRUE];if(ss=='pam50_basal')a=a[basal==TRUE];if(ss=='tnbc_MFS'){m=fread('results/tables/GSE58812_TNBC_survival_metadata.csv');idx=match(a$id,m$geo_accession);a[,time:=as.numeric(m$mfs_days[idx])];a[,event:=as.integer(m$meta[idx])]}
  a=a[complete.cases(a[,.(age,time,event,original,immune)]) & time>0]
  for(ad in c('age','immune')){
   fit=coxph(as.formula(paste('Surv(time,event)~original+age',if(ad=='immune')'+immune' else '')),data=a)
   hr=exp(coef(fit)['original']);ref=expected[dataset==ds & subset==ss & adjustment==ad]
   stopifnot(nrow(ref)==1,nrow(a)==ref$n,sum(a$event)==ref$events,abs(hr-ref$hr)<1e-8)
   checks[[length(checks)+1]]=data.table(dataset=ds,subset=ss,adjustment=ad,n=nrow(a),events=sum(a$event),hr_independent=hr,absolute_error=abs(hr-ref$hr))
  }
 }
}
fwrite(rbindlist(checks),file.path(t,'independent_cox_verification.csv'));writeLines(capture.output(sessionInfo()),'logs/revision_R1/verification_sessionInfo.txt');cat('Independent CSV-to-Cox checks passed:',length(checks),'\n')
