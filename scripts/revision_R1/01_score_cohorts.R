set.seed(123)
.libPaths(c('.Rlib',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(GSEABase);library(estimate)})
source('data/revision_R1/MCPcounter.R')
source('scripts/revision_R1/estimate_fast_io.R')
source('data/revision_R1/singscore_source/rankAndScoring.R')
source('data/revision_R1/singscore_source/simpleScoreGeneric.R')
T <- 'results/revision_R1/tables'; D <- 'data/revision_R1'
sets <- list(CD8=c('CD8A','CD8B','PDCD1','LAG3','HAVCR2','TIGIT','TOX','CXCL13','GZMB','PRF1','NKG7'),AP=c('HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','TAPBP','PSMB8','PSMB9','NLRC5'),IFN=c('IFNG','STAT1','IRF1','CXCL9','CXCL10','GBP1','GBP5','ISG15','IFIT1','MX1'))
writeLines(capture.output(sessionInfo()),'logs/revision_R1/scoring_sessionInfo.txt')
z <- function(x) {s=sd(x,na.rm=TRUE);if(!is.finite(s)||s==0)stop('Invalid score SD');(x-mean(x,na.rm=TRUE))/s}
mg <- as.data.frame(fread(file.path(D,'MCP_genes.txt'))); names(mg)<-c('HUGO symbols','Cell population','ENTREZID','ENSEMBL ID')
for(ds in c('GSE58812','GSE96058','METABRIC','TCGA')) {
 if(file.exists(file.path(D,paste0(ds,'_analysis.rds'))))next
 cat('START',ds,format(Sys.time()),'\n');flush.console()
 if(ds=='GSE58812') {
  e=fread('data/processed/GSE58812_GPL570_gene_symbol_expression.tsv.gz');c=fread('results/tables/GSE58812_TNBC_survival_metadata.csv')
  c[,`:=`(id=geo_accession,age=as.numeric(age_at_diag),time=as.numeric(os_days),event=as.integer(death),tnbc=TRUE,basal=NA,subtype='TNBC',stage=NA_character_)]; sam=c$id
 } else if(ds=='GSE96058') {
  e=fread('data/raw/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz');c=fread('results/tables/GSE96058_signature_survival_analysis_dataset.csv')
  c[,`:=`(id=primary_sample_id,age=as.numeric(age_at_diagnosis),time=as.numeric(os_days),event=as.integer(os_event),tnbc=fifelse(pathology_complete,tnbc_pathology,NA),basal=basal_pam50,subtype=pam50_subtype,stage=NA_character_)];sam=c$id
 } else {
  study=if(ds=='METABRIC')'brca_metabric' else 'brca_tcga_pan_can_atlas_2018'
  expr=if(ds=='METABRIC')'data_mrna_illumina_microarray.txt' else 'data_mrna_seq_v2_rsem.txt'
  e=fread(file.path(D,study,expr)); e[,Entrez_Gene_Id:=NULL]
  p=fread(file.path(D,study,'data_clinical_patient.txt'),skip=4,na.strings=c('NA','','[Not Available]'))
  s=fread(file.path(D,study,'data_clinical_sample.txt'),skip=4,na.strings=c('NA','','[Not Available]'))
  c=merge(s,p,by='PATIENT_ID');c=c[SAMPLE_TYPE=='Primary'];setorder(c,SAMPLE_ID);c=c[!duplicated(PATIENT_ID)]
  c[,`:=`(id=SAMPLE_ID,time=as.numeric(OS_MONTHS)*365.25/12,event=fifelse(OS_STATUS=='1:DECEASED',1L,fifelse(OS_STATUS=='0:LIVING',0L,NA_integer_)))]
  if(ds=='METABRIC') c[,`:=`(age=as.numeric(AGE_AT_DIAGNOSIS),tnbc=fifelse(ER_STATUS %in% c('Positive','Negative') & PR_STATUS %in% c('Positive','Negative') & HER2_STATUS %in% c('Positive','Negative'),ER_STATUS=='Negative' & PR_STATUS=='Negative' & HER2_STATUS=='Negative',NA),basal=fifelse(!is.na(CLAUDIN_SUBTYPE),CLAUDIN_SUBTYPE=='Basal',NA),subtype=CLAUDIN_SUBTYPE,stage=as.character(TUMOR_STAGE))]
  else c[,`:=`(age=as.numeric(AGE),tnbc=NA,basal=fifelse(!is.na(SUBTYPE),SUBTYPE=='BRCA_Basal',NA),subtype=sub('BRCA_','',SUBTYPE),stage=as.character(AJCC_PATHOLOGIC_TUMOR_STAGE))]
  sam=intersect(c$id,names(e)); c=c[id %in% sam]
 }
 setnames(e,1,'gene'); e[,gene:=toupper(gene)];e=e[!is.na(gene)&nzchar(gene)]
 sam=intersect(sam,names(e));stopifnot(!anyDuplicated(c$id),length(sam)>0)
 mat=as.matrix(e[,..sam]); storage.mode(mat)='double';rownames(mat)=e$gene
 if(anyDuplicated(e$gene)){valid=is.finite(mat); mm=mat;mm[!valid]=0;sm=rowsum(mm,e$gene,reorder=FALSE);ct=rowsum(valid*1,e$gene,reorder=FALSE);mat=sm/ct;rm(mm,sm,ct,valid)}
 rm(e);gc()
 if(ds=='TCGA')mat=log2(mat+1)
 # Do not impute matrix missingness. Only genes observed across the cohort enter the common ranking background.
 good=rowSums(!is.finite(mat))==0;removed=sum(!good);mat=mat[good,,drop=FALSE]
 cov=rbindlist(lapply(names(sets),function(k)data.table(dataset=ds,component=k,requested=length(sets[[k]]),observed=sum(sets[[k]]%in%rownames(mat)),missing=paste(setdiff(sets[[k]],rownames(mat)),collapse=';'))))
 fwrite(cov,file.path(T,paste0(ds,'_coverage.csv')));if(any(cov$observed!=cov$requested))stop(ds,' incomplete fixed signature')
 sg=mat[unlist(sets),,drop=FALSE];zg=t(apply(sg,1,z));components=sapply(sets,function(g)colMeans(zg[g,,drop=FALSE]))
 sc=data.table(id=colnames(mat),original=z(rowMeans(components)),gene_equal=z(colMeans(zg)))
 ranks=rankExpr(mat,tiesMethod='min');rs=sapply(sets,function(g)simpleScore(ranks,upSet=g)$TotalScore)
 sc[,singscore:=z(rowMeans(rs))];rm(ranks);gc()
 mc=t(MCPcounter.estimate(mat,'HUGO_symbols',genes=mg));sc=cbind(sc,as.data.table(mc))
 chosen=c('T cells','CD8 T cells','Cytotoxic lymphocytes','Monocytic lineage')
 overlap=rbindlist(lapply(chosen,function(k){g=mg[mg[['Cell population']]==k,'HUGO symbols'];data.table(dataset=ds,population=k,marker=g,present=g%in%rownames(mat),overlap=g%in%unlist(sets))}))
 fwrite(overlap,file.path(T,paste0(ds,'_MCP_marker_audit.csv')))
 for(k in chosen){g=setdiff(intersect(mg[mg[['Cell population']]==k,'HUGO symbols'],rownames(mat)),unlist(sets));sc[,(paste0(k,'_no_overlap')):=if(length(g))colMeans(mat[g,,drop=FALSE]) else NA_real_]}
 if(ds %in% c('GSE58812','GSE96058')) {
  est=fread('results/tables/estimate_scores.csv')[dataset==ds,.(id=sample_id,immune=ImmuneScore_z,stromal=StromalScore_z,purity=TumorPurity_z)]
  old=fread(paste0('results/tables/',ds,'_immune_reactivation_scores.csv'))
  oldid=if('sample_id'%in%names(old))'sample_id' else 'sample_title'
  check=merge(sc[,.(id,original)],old[,.(id=get(oldid),oldscore=immune_reactivation_score)],by='id');check[,oldscore:=z(oldscore)]
  fwrite(data.table(dataset=ds,max_abs_error=max(abs(check$original-check$oldscore))),file.path(T,paste0(ds,'_original_reproduction.csv')))
 } else {
  dir.create(file.path(D,'estimate'),showWarnings=FALSE)
  txt=file.path(D,'estimate',paste0(ds,'.txt'));gct=sub('.txt$','.gct',txt);out=sub('.txt$','_scores.gct',txt)
  if(!file.exists(out)){
   estimate_in_memory(mat,out,platform='illumina')
  }
  et=fread(out,skip=2);ev=as.matrix(et[,-c(1,2)]);rownames(ev)=et[[1]]
  est=data.table(id=sam[match(colnames(ev),make.names(sam))],immune=z(as.numeric(ev['ImmuneScore',])),stromal=z(as.numeric(ev['StromalScore',])),purity=NA_real_)
 }
 dat=merge(c,merge(sc,est,by='id'),by='id');dat[,dataset:=ds];stopifnot(nrow(dat)==nrow(c),!anyDuplicated(dat$id));saveRDS(dat,file.path(D,paste0(ds,'_analysis.rds')))
 fwrite(dat[,c('id','dataset','age','time','event','tnbc','basal','subtype','stage','original','gene_equal','singscore','immune','stromal','purity',chosen,paste0(chosen,'_no_overlap')),with=FALSE],file.path(T,paste0(ds,'_analysis_rows.csv')))
 fwrite(data.table(dataset=ds,expression_samples=ncol(mat),expression_genes=nrow(mat),nonfinite_genes_removed=removed,matched_samples=nrow(dat),os_complete=sum(complete.cases(dat[,.(time,event,age)]))),file.path(T,paste0(ds,'_input_qc.csv')))
 rm(mat,sg,zg);gc();cat('DONE',ds,format(Sys.time()),'\n');flush.console()
}
writeLines(capture.output(sessionInfo()),'logs/revision_R1/scoring_sessionInfo.txt')
