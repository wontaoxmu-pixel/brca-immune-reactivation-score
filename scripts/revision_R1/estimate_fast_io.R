# Same installed ESTIMATE algorithm; in-memory input adapter avoids slow read.table on large TSV.
# No scoring equations or reference sets changed. Validate against existing package outputs.
estimate_in_memory <- function(mat,out,platform='illumina'){
 common=get('common_genes',asNamespace('estimate'))
 genes=sort(intersect(common$GeneSymbol,rownames(mat)))
 df=data.frame(Description=genes,mat[genes,,drop=FALSE],check.names=TRUE,row.names=genes)
 fn=estimate::estimateScore
 en=new.env(parent=environment(fn));en$read.delim=function(...)df
 environment(fn)=en
 fn('in_memory_common_genes.gct',out,platform=platform)
 invisible(length(genes))
}
