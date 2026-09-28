set.seed(123);.libPaths(c('.Rlib',.libPaths()));suppressPackageStartupMessages({library(data.table);library(survival)})
sets=list(c('GSE58812','TNBC'),c('GSE96058','All primary'),c('GSE96058','Pathology TNBC'),c('GSE96058','Basal-like'));rows=list()
draw=function(){par(mfrow=c(2,2),mar=c(4,4,2.5,1),mgp=c(2.2,.7,0),cex=.85)
 for(i in seq_along(sets)){
  ds=sets[[i]][1];ss=sets[[i]][2];d=readRDS(paste0('data/revision_R1/',ds,'_analysis.rds'))
  if(ss=='Pathology TNBC')d=d[tnbc==TRUE];if(ss=='Basal-like')d=d[basal==TRUE]
  d=d[complete.cases(d[,.(time,event,original)]) & time>0];d[,group:=factor(ifelse(original>median(original),'Above median','At/below median'),levels=c('At/below median','Above median'))]
  fit=survfit(Surv(time/365.25,event)~group,data=d)
  plot(fit,col=c('#C67B21','#2878B5'),lwd=1.5,conf.int=FALSE,mark.time=TRUE,xlab='Time (years)',ylab='Overall survival probability',main=paste0(LETTERS[i],'. ',ds,' / ',ss),ylim=c(0,1));legend('bottomleft',legend=paste0(levels(d$group),' (n=',as.integer(table(d$group)),')'),col=c('#C67B21','#2878B5'),lty=1,bty='n',cex=.75)
  s=summary(fit);rows[[i]]<<-data.table(dataset=ds,subset=ss,time_years=s$time,survival=s$surv,lower=s$lower,upper=s$upper,n_risk=s$n.risk,n_event=s$n.event,stratum=as.character(s$strata))
 }
}
png('results/revision_R1/figures/S4_Fig.png',width=2100,height=1600,res=300);draw();dev.off()
pdf('results/revision_R1/figures/S4_Fig.pdf',width=7,height=5.33);draw();dev.off()
tiff('results/revision_R1/figures/S4_Fig.tif',width=4200,height=3200,res=600,compression='lzw');draw();dev.off()
fwrite(rbindlist(rows),'results/revision_R1/tables/exploratory_KM_source.csv');writeLines(capture.output(sessionInfo()),'logs/revision_R1/KM_sessionInfo.txt')
