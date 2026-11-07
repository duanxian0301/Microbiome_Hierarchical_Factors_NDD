if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
root <- file.path(Sys.getenv("MHFNDD_PROJECT_ROOT"), "FULL_SCALE_GFA_FINAL_EVALUATION")
out <- file.path(root,'models'); dir.create(out,showWarnings=FALSE,recursive=TRUE)
meta <- jsonlite::fromJSON(file.path(root,'r241_inputs','R241_input_metadata.json')); K<-meta$K_traits; traits<-as.character(meta$trait_order); training_J<-330704L
Z<-matrix(readBin(file.path(root,'HUNT241_genomewide_ld_controlled_Z.float64.bin'),what='double',n=training_J*K,size=8,endian='little'),nrow=training_J,ncol=K)
if(length(Z)!=training_J*K||any(!is.finite(Z)))stop('Invalid frozen HUNT training Z matrix'); colnames(Z)<-traits
R<-as.matrix(read.delim(file.path(root,'HUNT_R241_ldsc_raw_matrix.tsv'),check.names=FALSE));storage.mode(R)<-'double';rownames(R)<-colnames(R)<-traits
if(!identical(dim(R),c(K,K))||any(!is.finite(R))||max(abs(R-t(R)))>1e-12||min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)<=0)stop('Invalid HUNT raw R_ldsc input')
set.seed(20261014);started<-Sys.time();fit<-GFA::gfa_fit(Z_hat=Z,N=setNames(rep(12652,K),traits),R=R);ended<-Sys.time()
saveRDS(fit,file.path(out,'HUNT241_genomewide_ldcontrolled_Rldsc_seed20261014.rds'))
write.table(data.frame(started=as.character(started),ended=as.character(ended),elapsed_seconds=as.numeric(difftime(ended,started,units='secs')),training_snps=training_J,traits=K,retained_data_factors=ncol(fit$F_hat)),file.path(out,'HUNT241_genomewide_ldcontrolled_Rldsc_seed20261014_timing.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'HUNT241_genomewide_ldcontrolled_Rldsc_seed20261014_sessionInfo.txt'))
writeLines(c(sprintf('GFA_path=%s',find.package('GFA')),sprintf('GFA_version=%s',as.character(packageVersion('GFA'))),sprintf('flashier_version=%s',as.character(packageVersion('flashier')))),file.path(out,'HUNT241_genomewide_ldcontrolled_Rldsc_seed20261014_software.txt'))
