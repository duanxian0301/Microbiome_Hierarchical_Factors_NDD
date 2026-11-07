if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
root <- file.path(Sys.getenv("MHFNDD_PROJECT_ROOT"), "PHASE14_TARGETED_VALIDATION", "B_SCZ_CD_UC_LOCI")
indir <- file.path(root,"coloc_inputs")
files <- list.files(indir,pattern="__coloc_input\\.tsv$",full.names=TRUE)
if (!length(files)) stop("No coloc inputs")
rows <- list()
for (f in files) {
  x <- read.delim(f,check.names=FALSE,stringsAsFactors=FALSE)
  x <- x[is.finite(x$DISEASE_BETA) & is.finite(x$DISEASE_SE) & x$DISEASE_SE>0 &
         is.finite(x$SPECIES_BETA_ALIGNED) & is.finite(x$SPECIES_SE) & x$SPECIES_SE>0 &
         is.finite(x$REFERENCE_MAF_1000G_EUR) & x$REFERENCE_MAF_1000G_EUR>0 & x$REFERENCE_MAF_1000G_EUR<0.5,]
  region <- sub("__.*$","",basename(f))
  disease <- sub("_.*$","",region)
  sp <- if (nrow(x)) x$SPECIES[[1]] else sub("^[^_]+__|__coloc_input.tsv$","",basename(f))
  if (nrow(x)<100) {
    rows[[length(rows)+1]] <- data.frame(REGION_ID=region,DISEASE=disease,SPECIES=sp,
      COMMON_SNP_N=nrow(x),METHOD="coloc.abf",STATUS="NOT_RUN_INSUFFICIENT_COMMON_SNPS",
      PP_H0=NA,PP_H1=NA,PP_H2=NA,PP_H3=NA,PP_H4=NA,PP_H4_P12_1E6=NA,PP_H4_P12_1E4=NA,COLOC_TOP_SNP=NA,COLOC_TOP_SNP_PP_H4=NA,
      SPECIES_MIN_P=if(nrow(x)) min(x$SPECIES_P,na.rm=TRUE) else NA,
      DISEASE_MIN_P=if(nrow(x)) min(x$DISEASE_P,na.rm=TRUE) else NA,
      INPUT_LIMITATION="fewer than 100 harmonized common variants",stringsAsFactors=FALSE)
    next
  }
  ncase <- median(x$DISEASE_NCASE,na.rm=TRUE); ncontrol <- median(x$DISEASE_NCONTROL,na.rm=TRUE)
  d1 <- list(beta=x$SPECIES_BETA_ALIGNED,varbeta=x$SPECIES_SE^2,snp=x$SNP_KEY,
             position=x$BP,MAF=x$REFERENCE_MAF_1000G_EUR,N=16017,type="quant")
  d2 <- list(beta=x$DISEASE_BETA,varbeta=x$DISEASE_SE^2,snp=x$SNP_KEY,
             position=x$BP,MAF=x$REFERENCE_MAF_1000G_EUR,N=ncase+ncontrol,
             s=ncase/(ncase+ncontrol),type="cc")
  ans <- tryCatch(coloc::coloc.abf(d1,d2),error=function(e)e)
  if (inherits(ans,"error")) {
    rows[[length(rows)+1]] <- data.frame(REGION_ID=region,DISEASE=disease,SPECIES=sp,
      COMMON_SNP_N=nrow(x),METHOD="coloc.abf",STATUS=paste0("ERROR: ",conditionMessage(ans)),
      PP_H0=NA,PP_H1=NA,PP_H2=NA,PP_H3=NA,PP_H4=NA,PP_H4_P12_1E6=NA,PP_H4_P12_1E4=NA,COLOC_TOP_SNP=NA,COLOC_TOP_SNP_PP_H4=NA,
      SPECIES_MIN_P=min(x$SPECIES_P),DISEASE_MIN_P=min(x$DISEASE_P),
      INPUT_LIMITATION="reference MAF and study-level species N sensitivity; original species INFO/SNP-N absent",stringsAsFactors=FALSE)
  } else {
    s <- ans$summary
    low <- coloc::coloc.abf(d1,d2,p12=1e-6)$summary
    high <- coloc::coloc.abf(d1,d2,p12=1e-4)$summary
    topi <- which.max(ans$results$SNP.PP.H4)
    rows[[length(rows)+1]] <- data.frame(REGION_ID=region,DISEASE=disease,SPECIES=sp,
      COMMON_SNP_N=as.integer(s[["nsnps"]]),METHOD="coloc.abf",STATUS="COMPLETED_SENSITIVITY",
      PP_H0=as.numeric(s[["PP.H0.abf"]]),PP_H1=as.numeric(s[["PP.H1.abf"]]),
      PP_H2=as.numeric(s[["PP.H2.abf"]]),PP_H3=as.numeric(s[["PP.H3.abf"]]),
      PP_H4=as.numeric(s[["PP.H4.abf"]]),PP_H4_P12_1E6=as.numeric(low[["PP.H4.abf"]]),PP_H4_P12_1E4=as.numeric(high[["PP.H4.abf"]]),COLOC_TOP_SNP=as.character(ans$results$snp[[topi]]),COLOC_TOP_SNP_PP_H4=as.numeric(ans$results$SNP.PP.H4[[topi]]),SPECIES_MIN_P=min(x$SPECIES_P),DISEASE_MIN_P=min(x$DISEASE_P),
      INPUT_LIMITATION="1000G EUR reference MAF; Swedish study-level N=16017; original species INFO and SNP-specific N unavailable; disease INFO unavailable for CD/UC",
      stringsAsFactors=FALSE)
  }
}
out <- do.call(rbind,rows)
write.table(out,file.path(root,"PHASE14B_COLOCALIZATION_RESULTS.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
jsonlite::write_json(list(status="complete",models=nrow(out),completed=sum(out$STATUS=="COMPLETED_SENSITIVITY")),file.path(root,"checkpoints","coloc_checkpoint.json"),pretty=TRUE,auto_unbox=TRUE)
print(out)
