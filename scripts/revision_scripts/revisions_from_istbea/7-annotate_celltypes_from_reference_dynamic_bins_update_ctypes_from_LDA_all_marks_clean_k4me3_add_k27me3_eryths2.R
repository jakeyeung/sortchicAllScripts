# Jake Yeung
# Date of Creation: 2022-04-13
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/7-annotate_celltypes_from_reference_dynamic_bins_update_ctypes_from_LDA_all_marks_clean_k4me3_add_k27me3_eryths.R
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

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
# jmarks <- c("k27me3"); names(jmarks) <- jmarks
jmark <- jmarks[[1]]

jmarks.active <- c("k4me1", "k4me3"); names(jmarks.active)
jmarks.repress <- c("k27me3", "k9me3"); names(jmarks.repress)

jsuffix <- "dynamicbins"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"
dir.create(outdir)
# for (jmark in jmarks){

ctypes.to.annotate <- c("Basophils", "Bcells", "CMP", "DCs", "Eryths", "Granulocytes", "HSCs", "MEP", "Monocytes", "MPPs", "NKs", "pDCs", "ST", "LT")
names(ctypes.to.annotate) <- ctypes.to.annotate

LL.all.lst <- mclapply(jmarks, function(jmark){
  print(jmark)
  
  outrds <- file.path(outdir, paste0("LLmat_by_batch_", jsuffix, "_", jmark, ".", Sys.Date(), ".rds"))
  
  # Load metadata -----------------------------------------------------------
  
  if (jmark == "k4me1"){
    # inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/varfilt/metadata_celltyping_k4me1.allbins.2022-02-01.txt"
    inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_with_colors_keepTcells.k4me1.2022-04-13.txt"
    # indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
    # inf.lda <- ""
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
    # inf.lda <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata_k4me3_cleaned_up/metadata_k4me3_louvain_filtered.no3no7.2022-04-12.txt"
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_clean_eryths2"
    inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".dynamicbins.2022-04-16.txt"))
    # inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-15/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned/metadata_plate_experi_batch.cleaned.dynamicbins.k9me3.2022-02-17.txt"
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  
  load(inf.ldaout, v=T)
  dat.meta <- fread(inf.meta)
  print(unique(sort(dat.meta$ctype)))
  assertthat::assert_that(nrow(dat.meta) == ncol(count.mat))
  
  ggplot(dat.meta, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # reclassify celltypes for AllCells, HSPCs, IL7RLinNeg, LinNeg, LSK
  
  # Load rawcounts ----------------------------------------------------------

  
  assertthat::assert_that(ncol(count.mat) == nrow(dat.meta))
  tm <- posterior(out.lda)
  dat.impute.linear <- t(tm$topics %*% tm$terms)
  
  # Annotate celltypes ------------------------------------------------------
  
  dat.meta.ctype.filt <- subset(dat.meta, ctype %in% ctypes.to.annotate)
  dat.meta.testset <- subset(dat.meta, !ctype %in% ctypes.to.annotate)
  cluster.ids <- lapply(split(dat.meta.ctype.filt, f = dat.meta.ctype.filt$ctype), function(x) x$cell)
  
  # for DCs: use only new sorts
  cluster.ids.orig <- cluster.ids
  #print(lapply(cluster.ids.orig, length))
  cluster.ids$DCs <- cluster.ids$DCs[grepl("^PZ-sortChIC-BM", cluster.ids$DCs)]  # 617 -> 125 cells
  # use Bcells
  cluster.ids$Bcells <- cluster.ids$Bcells[grepl("^PZ-sortChIC-BM", cluster.ids$Bcells)]
  cluster.ids$Bcells.old <- cluster.ids.orig$Bcells[!grepl("^PZ-sortChIC-BM", cluster.ids.orig$Bcells)]
  cluster.ids$Eryths <- cluster.ids$Eryths[grepl("^PZ-sortChIC-BM", cluster.ids$Eryths)]
  cluster.ids$Eryths.old <- cluster.ids.orig$Eryths[!grepl("^PZ-sortChIC-BM", cluster.ids.orig$Eryths)]
  cluster.ids$Granulocytes <- cluster.ids$Granulocytes[grepl("^PZ-sortChIC-BM", cluster.ids$Granulocytes)]
  cluster.ids$Granulocytes.old <- cluster.ids.orig$Granulocytes[!grepl("^PZ-sortChIC-BM", cluster.ids.orig$Granulocytes)]
  
  # for NKs: use old data
  # cluster.ids$NKs <- cluster.ids$NKs[grepl("^PZ-ChIC", cluster.ids$NKs)]  # 617 -> 125 cells
  if (jmark != "k9me3"){
    # bug otherwise because k9me3 old has no NK cells
    cluster.ids$NKs <- cluster.ids$NKs[!grepl("^PZ-sortChIC-BM", cluster.ids$NKs)]  # 435 -> 160 cells
  }
  # jtest <- cluster.ids$NKs[!grepl("^PZ-sortChIC-BM", cluster.ids$NKs)]  # 435 -> 160 cells
  # for 
  # lapply(cluster.ids, length)
  
  
  
  # count.mat.pbulk <- SumAcrossClusters(count.mat = count.mat, cnames.keep.lst = cluster.ids)
  # count.mat.pbulk <- do.call(cbind, count.mat.pbulk)
  # 
  # #' zero counts are problematic when calculating likelihoods, add pseudocount to compensate
  # mats.keep.pseudobulk.pcount <- count.mat.pbulk + 1
  # mats.keep.pseudobulk.pcount.frac <- sweep(mats.keep.pseudobulk.pcount, MARGIN = 2, STATS = colSums(mats.keep.pseudobulk.pcount), FUN = "/")
  
  ref.probs.mat <- as.data.frame(ApplyAcrossClusters(count.mat = dat.impute.linear, cnames.keep.lst = cluster.ids, fn = mean))
  
  
  # Estimate from multionomials  --------------------------------------------
  
  # (cellname <- dat.meta.testset$cell[[1]])
  # jcell <- count.mat[, cellname]
  # 
  # #' note the input here is unnormalized counts. No need for log-transform or scaling rows or columns.
  # #' should be robust to downsampling too.
  # LL <- AnnotateCelltypes::CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
  
  # run on all
  LL.all <- apply(count.mat, 2, function(jcell){
    LL <- AnnotateCelltypes::CalculateMultinomLL(jcell, ref.probs.mat = ref.probs.mat)
  })
  rownames(LL.all) <- colnames(ref.probs.mat)
  saveRDS(LL.all, file = outrds)
  print(paste("Done for", jmark))
  return(LL.all)
  # }
}, mc.cores = length(jmarks))
  
  
  


