if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop('usage: cohort chunk_index chunks_total')
cohort <- args[[1]]; chunk_i <- as.integer(args[[2]]); chunks <- as.integer(args[[3]])
root <- Sys.getenv("MHFNDD_PROJECT_ROOT")
out <- file.path(root,'FULL_SCALE_GFA_FINAL_EVALUATION'); inp <- file.path(out,'r241_inputs')
gfa <- file.path(root,'GFA_NATIVE_MVMR_TRANSFER_AUDIT','official_GFA','R')
suppressPackageStartupMessages({library(dplyr); library(purrr); library(Matrix)})
source(file.path(gfa,'ldsc.R')); source(file.path(gfa,'ldsc_rg.R')); source(file.path(gfa,'estimate_R.R'))
meta <- jsonlite::fromJSON(file.path(inp,'R241_input_metadata.json'))
J <- meta$J_common_allele_verified; traits <- as.character(meta$trait_order); K <- length(traits)
zpath <- file.path(inp,paste0(cohort,'_Z_ldsc.float64.bin'))
started <- Sys.time()
Z <- matrix(readBin(zpath,what='double',n=J*K,size=8,endian='little'),nrow=J,ncol=K)
if(length(Z)!=J*K || any(!is.finite(Z))) stop('invalid Z input')
colnames(Z) <- traits
ld <- readBin(file.path(inp,'hm3_nomHC_L2.float64.bin'),what='double',n=J,size=8,endian='little')
if(length(ld)!=J || any(!is.finite(ld)) || any(ld<=0)) stop('invalid LD score input')
N0 <- if(cohort=='Swedish') 16017 else if(cohort=='HUNT') 12652 else stop('unknown cohort')
N <- setNames(rep(N0,K),traits)
pairs <- data.frame(trait1=rep(traits,times=seq_len(K)),trait2=traits[sequence(seq_len(K))])
take <- ((seq_len(nrow(pairs))-1) %% chunks)+1L == chunk_i
pairs <- pairs[take,,drop=FALSE]
ans <- R_ldsc(Z_hat=Z,ldscores=ld,ld_size=meta$ld_size_sum_M_official_weight_reference,
              N=N,return_gencov=FALSE,make_well_conditioned=FALSE,blocks=NULL,ncores=1,
              comparisons=pairs)
ans$cohort <- cohort; ans$common_allele_verified_snps <- J; ans$ld_size_M <- meta$ld_size_sum_M_official_weight_reference
ans$N_trait1 <- N0; ans$N_trait2 <- N0; ans$official_method <- 'GFA::R_ldsc raw intercept, make_well_conditioned=FALSE, blocks=NULL'
ans$chunk_index <- chunk_i; ans$chunks_total <- chunks
raw <- file.path(out,'r241_raw'); dir.create(raw,showWarnings=FALSE,recursive=TRUE)
write.table(ans,file.path(raw,sprintf('%s_R_ldsc_chunk_%02d_of_%02d.tsv',cohort,chunk_i,chunks)),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(raw,sprintf('%s_R_ldsc_chunk_%02d_sessionInfo.txt',cohort,chunk_i)))
writeLines(sprintf('started=%s\nended=%s\nelapsed_seconds=%.3f\npairs=%d\n',started,Sys.time(),as.numeric(difftime(Sys.time(),started,units='secs')),nrow(pairs)),file.path(raw,sprintf('%s_R_ldsc_chunk_%02d_timing.txt',cohort,chunk_i)))
