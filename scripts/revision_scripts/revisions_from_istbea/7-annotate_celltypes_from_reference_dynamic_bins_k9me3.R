# Jake Yeung
# Date of Creation: 2022-02-01
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/7-annotate_celltypes_from_reference.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)
library(AnnotateCelltypes)
library(parallel)

jmarks <- c("k9me3"); names(jmarks) <- jmarks
jmark <- jmarks[[1]]

jsuffix <- "dynamicbins"
# outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_newonly"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping"
dir.create(outdir)
# for (jmark in jmarks){

LL.all.lst <- mclapply(jmarks, function(jmark){
  print(jmark)
  
  outrds <- file.path(outdir, paste0("LLmat_", jsuffix, "_", jmark, ".", Sys.Date(), ".rds"))
  
  # Load metadata -----------------------------------------------------------
  
  indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/varfilt"
  inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".allbins.2022-02-01.txt"))
  dat.meta <- fread(inf.meta)
  
  print(unique(sort(dat.meta$ctype)))
  
  # reclassify celltypes for AllCells, HSPCs, IL7RLinNeg, LinNeg, LSK
  
  
  # Load rawcounts ----------------------------------------------------------
  
  indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
  inf.ldaout <- file.path(indir.ldaout, paste0("ldaOut.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj"))
  load(inf.ldaout, v=T)
  
  
  
  # Annotate celltypes ------------------------------------------------------
  
  # ctypes.to.annotate <- c("Basophils", "Bcells", "CMP", "DCs", "Eryths", "GMP", "Granulocytes", "HSCs", "LT", "MEP", "Monocytes", "MPPs", "NKs", "pDCs", "ST", "Tcells")
  ctypes.to.annotate <- c("Bcells", "Eryths", "GMP", "DCs", "Granulocytes", "HSCs", "LT", "MEP", "MPPs", "pDCs", "ST", "Monocytes")
  names(ctypes.to.annotate) <- ctypes.to.annotate
  
  # dat.meta.ctype.filt <- subset(dat.meta, ctype %in% ctypes.to.annotate & batch == "New")
  dat.meta.ctype.filt <- subset(dat.meta, ctype %in% ctypes.to.annotate)
  dat.meta.testset <- subset(dat.meta, !ctype %in% ctypes.to.annotate)
  cluster.ids <- lapply(split(dat.meta.ctype.filt, f = dat.meta.ctype.filt$ctype), function(x) x$cell)
  
  
  count.mat.pbulk <- SumAcrossClusters(count.mat = count.mat, cnames.keep.lst = cluster.ids)
  count.mat.pbulk <- do.call(cbind, count.mat.pbulk)
  
  
  #' zero counts are problematic when calculating likelihoods, add pseudocount to compensate
  mats.keep.pseudobulk.pcount <- count.mat.pbulk + 1
  mats.keep.pseudobulk.pcount.frac <- sweep(mats.keep.pseudobulk.pcount, MARGIN = 2, STATS = colSums(mats.keep.pseudobulk.pcount), FUN = "/")
  
  # Estimate from multionomials  --------------------------------------------
  
  # (cellname <- dat.meta.testset$cell[[1]])
  # jcell <- count.mat[, cellname]
  # 
  # #' note the input here is unnormalized counts. No need for log-transform or scaling rows or columns.
  # #' should be robust to downsampling too.
  # LL <- AnnotateCelltypes::CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
  
  # run on all
  LL.all <- apply(count.mat, 2, function(jcell){
    LL <- AnnotateCelltypes::CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
  })
  rownames(LL.all) <- colnames(mats.keep.pseudobulk.pcount.frac)
  saveRDS(LL.all, file = outrds)
  print(paste("Done for", jmark))
  return(LL.all)
  # }
}, mc.cores = length(jmarks))
  
  
  


