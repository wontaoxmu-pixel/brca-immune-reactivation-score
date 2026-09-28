set.seed(123);.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(survival)})
T='results/revision_R1/tables';D='data/revision_R1';B=2000L; setDTthreads(1)
z=function(x)as.numeric(scale(x)); rows=list(); diagrows=list(); incr=list(); inter=list(); precision=list(); bootout=list(); failures=list(); corr=list()
fit=function(f,d)tryCatch(coxph(f,data=d,x=TRUE,y=TRUE,model=TRUE,control=coxph.control(iter.max=50)),warning=function(w)NULL,error=function(e)NULL)
ci=function(m)unname(summary(m)$concordance[1])
term=function(m,k){if(is.null(m))return(NULL);s=summary(m);if(!k%in%rownames(s$coefficients))return(NULL);data.table(hr=s$conf.int[k,1],lower=s$conf.int[k,3],upper=s$conf.int[k,4],p=s$coefficients[k,5])}
for(ds in c('GSE58812','GSE96058','METABRIC','TCGA')){
 ready=file.path(D,paste0(ds,'_analysis.rds')); waited=0;while(!file.exists(ready)&&waited<1800){Sys.sleep(5);waited=waited+5};if(!file.exists(ready))stop('Scoring not ready: ',ds)
 all=readRDS(file.path(D,paste0(ds,'_analysis.rds')));if(ds%in%c('METABRIC','TCGA'))all[,purity:=NA_real_];cat('START',ds,'\n');flush.console()
 abund=c('immune','T cells','CD8 T cells','Cytotoxic lymphocytes','Monocytic lineage')
 for(sc in c('original','gene_equal','singscore'))for(ab in c(abund,paste0(abund[-1],'_no_overlap'))){
  ok=complete.cases(all[,c(sc,ab),with=FALSE]);a=all[[sc]][ok];b=all[[ab]][ok];if(length(a)<20||sd(b)==0)next
  rr=replicate(B,{j=sample.int(length(a),replace=TRUE);cor(a[j],b[j],method='spearman')})
  corr[[length(corr)+1]]=data.table(dataset=ds,score=sc,abundance=ab,n=length(a),rho=cor(a,b,method='spearman'),lower=quantile(rr,.025,na.rm=TRUE),upper=quantile(rr,.975,na.rm=TRUE))
 }
 specs=list(all_primary=all);if(ds=='GSE58812')specs=list(tnbc=all)
 else{if(any(all$tnbc %in% TRUE))specs$pathology_tnbc=all[tnbc %in% TRUE];if(any(all$basal %in% TRUE))specs$pam50_basal=all[basal %in% TRUE]}
 if(ds=='GSE58812'){m=copy(all);m[,`:=`(time=as.numeric(mfs_days),event=as.integer(meta))];specs$tnbc_MFS=m}
 for(ss in names(specs)){
  d=copy(specs[[ss]]);d=d[complete.cases(d[,.(time,event,age,original,immune)]) & time>0 & event%in%c(0,1)]
  n=nrow(d);ev=sum(d$event);if(n<40||ev<20){failures[[length(failures)+1]]=data.table(dataset=ds,subset=ss,n=n,events=ev,reason='below n40 or 20 events');next}
  for(v in c(1,.5,.25))precision[[length(precision)+1]]=data.table(dataset=ds,subset=ss,n=n,events=ev,residual_variance=v,detectable_HR=exp((qnorm(.975)+qnorm(.8))/sqrt(ev*v)),protective_HR=exp(-(qnorm(.975)+qnorm(.8))/sqrt(ev*v)))
  for(sc in c('original','gene_equal','singscore')){
   d[,score:=get(sc)]
   for(adj in c('age','immune','T cells','CD8 T cells','Cytotoxic lymphocytes','Monocytic lineage','within_subset_SD','purity','stromal','stage')){
    q=copy(d);rhs='score+age'
    if(adj=='within_subset_SD')q[,score:=z(score)]
    if(adj%in%c(abund,'purity','stromal')){q[,ab:=z(get(adj))];rhs='score+age+ab'}
    if(adj=='stage'){
     q=q[!is.na(stage)&!stage%in%c('STAGE X','0','Not Available')];q[,stagef:=factor(stage)];if(nrow(q)<40||nlevels(q$stagef)<2)next
     rhs='score+age+stagef';if(sum(q$event)<10*(2+nlevels(q$stagef)-1))next
    }
    q=q[complete.cases(q[,unique(c('time','event','age','score',if(grepl('ab',rhs,fixed=TRUE))'ab')),with=FALSE])]
    m=fit(as.formula(paste('Surv(time,event)~',rhs)),q);r=term(m,'score')
    if(is.null(r)){failures[[length(failures)+1]]=data.table(dataset=ds,subset=ss,n=nrow(q),events=sum(q$event),reason=paste('fit failed',sc,adj));next}
    r[,`:=`(dataset=ds,subset=ss,score=sc,adjustment=adj,n=nrow(q),events=sum(q$event),family=if(sc=='original')'original_prognosis' else 'alternative_scores')];rows[[length(rows)+1]]=r
    ph=tryCatch(cox.zph(m)$table,error=function(e)NULL)
    if(!is.null(ph)){
     diagrows[[length(diagrows)+1]]=data.table(dataset=ds,subset=ss,score=sc,adjustment=adj,term=rownames(ph),ph_p=ph[,'p'])
     if(sc=='original' && adj=='age' && ph['score','p']<.05){
      mt=tryCatch(coxph(Surv(time,event)~score+age+tt(score),data=q,tt=function(x,t,...)x*log(pmax(t,1)/365.25)),error=function(e)NULL)
      rt=term(mt,'tt(score)');if(!is.null(rt)){rt[,`:=`(dataset=ds,subset=ss,score=sc,adjustment='score_by_log_time',n=nrow(q),events=sum(q$event),family='PH_sensitivity')];rows[[length(rows)+1]]=rt}
     }
    }
   }
   forms=list(age=Surv(time,event)~age,score=Surv(time,event)~age+score,immune=Surv(time,event)~age+immune,both=Surv(time,event)~age+score+immune)
   ms=lapply(forms,fit,d=d);if(any(vapply(ms,is.null,logical(1))))next
   ir=data.table(dataset=ds,subset=ss,score=sc,n=n,events=ev,c_age=ci(ms$age),c_score=ci(ms$score),c_immune=ci(ms$immune),c_both=ci(ms$both),score_added_p=pchisq(2*(ms$both$loglik[2]-ms$immune$loglik[2]),1,lower.tail=FALSE),immune_added_p=pchisq(2*(ms$both$loglik[2]-ms$score$loglik[2]),1,lower.tail=FALSE),AIC_score=AIC(ms$score),AIC_immune=AIC(ms$immune),AIC_both=AIC(ms$both));incr[[length(incr)+1]]=ir
   if(sc!='original')next
   cat('BOOT',ds,ss,n,ev,'\n');flush.console()
   # Pre-generate paired draws for deterministic parallelism. Every model uses the same draw.
   draws=replicate(B,sample.int(n,n,replace=TRUE),simplify=FALSE)
   one=function(j){dd=d[j];if(sum(dd$event)<5)return(rep(NA_real_,8));mm=lapply(forms,fit,d=dd);if(any(vapply(mm,is.null,logical(1))))return(rep(NA_real_,8));tr=vapply(mm,ci,numeric(1));te=vapply(mm,function(m){lp=as.numeric(predict(m,newdata=d,type='lp'));concordance(Surv(d$time,d$event)~lp,reverse=TRUE)$concordance},numeric(1));c(tr,tr-te)}
   vals=lapply(draws,function(j)tryCatch(one(j),error=function(e)rep(NA_real_,8)));v=do.call(rbind,vals);valid=complete.cases(v);v=v[valid,,drop=FALSE];nv=nrow(v)
   saveRDS(list(values=v,valid=valid,columns=c(paste0('apparent_',names(forms)),paste0('optimism_',names(forms)))),file.path(D,paste0('bootstrap_',ds,'_',ss,'.rds')))
   point=vapply(ms,ci,numeric(1))
   for(k in 1:4){qs=if(nv>=.9*B)quantile(v[,k],c(.025,.975)) else c(NA,NA);bootout[[length(bootout)+1]]=data.table(dataset=ds,subset=ss,metric=names(forms)[k],estimate=point[k],lower=qs[1],upper=qs[2],optimism=mean(v[,k+4]),corrected=point[k]-mean(v[,k+4]),valid=nv,attempted=B,n=n,events=ev,CI_target='refitted apparent C-index; corrected point separately')}
   dif=v[,4]-v[,3];op=v[,8]-v[,7];qs=if(nv>=.9*B)quantile(dif,c(.025,.975)) else c(NA,NA)
   bootout[[length(bootout)+1]]=data.table(dataset=ds,subset=ss,metric='delta_both_minus_immune',estimate=point[4]-point[3],lower=qs[1],upper=qs[2],optimism=mean(op),corrected=point[4]-point[3]-mean(op),valid=nv,attempted=B,n=n,events=ev,CI_target='paired refitted apparent delta; corrected point separately')
  }
 }
 for(groupvar in c('basal','tnbc')){
  d=all[complete.cases(all[,c('time','event','age','original',groupvar),with=FALSE])&time>0];d[,group:=as.integer(get(groupvar))];d[,score:=original]
  if(nrow(d)<40||length(unique(d$group))<2||any(tapply(d$event,d$group,sum)<20))next
  m0=fit(Surv(time,event)~score+age+group,d);m1=fit(Surv(time,event)~score*group+age,d)
  if(!is.null(m0)&&!is.null(m1))inter[[length(inter)+1]]=data.table(dataset=ds,group_definition=groupvar,n=nrow(d),events=sum(d$event),p=pchisq(2*(m1$loglik[2]-m0$loglik[2]),1,lower.tail=FALSE))
 }
 fwrite(rbindlist(corr,fill=TRUE),file.path(T,'abundance_correlations.csv'));fwrite(rbindlist(bootout,fill=TRUE),file.path(T,'bootstrap_performance.csv'))
 cat('DONE',ds,format(Sys.time()),'\n');flush.console()
}
a=rbindlist(rows,fill=TRUE);a[,q:=p.adjust(p,'BH'),by=family];fwrite(a,file.path(T,'cox_results.csv'))
b=rbindlist(incr,fill=TRUE);qq=p.adjust(c(b$score_added_p,b$immune_added_p),'BH');b[,score_added_q:=qq[seq_len(.N)]];b[,immune_added_q:=qq[nrow(b)+seq_len(.N)]];fwrite(b,file.path(T,'incremental_results.csv'))
c=rbindlist(inter,fill=TRUE);if(nrow(c))c[,q:=p.adjust(p,'BH')];fwrite(c,file.path(T,'interaction_results.csv'))
fwrite(rbindlist(diagrows,fill=TRUE),file.path(T,'PH_diagnostics.csv'));fwrite(rbindlist(precision,fill=TRUE),file.path(T,'precision_scenarios.csv'));fwrite(rbindlist(failures,fill=TRUE),file.path(T,'model_failures.csv'))
writeLines(capture.output(sessionInfo()),'logs/revision_R1/survival_sessionInfo.txt')
