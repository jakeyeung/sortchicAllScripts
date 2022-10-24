# Jake Yeung
# Date of Creation: 2022-05-02
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/23-correct_batch_effect_other_marks.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

jstart <- Sys.time() 

# jmark <- "k27me3"
# jmarkold <- "H3K27me3"

# jmarks <- c("k4me1", "k4me3", "k9me3", "k27me3"); names(jmarks) <- jmarks
jmarks <- c("k27me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


dat.meta.lst <- lapply(jmarks, function(jmark){
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta)
})


for (jmark in jmarks){
  
  print(jmark)
  
  outrds <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_", jmark, "_batch_corrected.", Sys.Date(), ".rds")
  outrds.lst <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_", jmark, "_batch_corrected_lst.", Sys.Date(), ".rds")
  outrdsraw <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_", jmark, "_raw_counts.", Sys.Date(), ".rds")
  
  if (jmark == "k4me1"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    
    jsuffix2 <- "merged_with_old_dynbins"
    indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_", jmark, "_clean_eryths")
    dname <- paste0("count_mat_", jsuffix2, ".", jmark, ".2022-04-15")
    inf.ldaout <- file.path(indir, paste0("lda_outputs.", dname), paste0("ldaOut.", dname, ".Robj"))
    
    # inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_new_only_dynbins.k27me3.2022-04-15/ldaOut.count_mat_new_only_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  
  load(inf.ldaout, v=T)
  tm <- posterior(out.lda)
  
  dat.impute.log2 <- t(log2(tm$topics %*% tm$terms))
  
  # Load meta ---------------------------------------------------------------
  
  # indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"
  # inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.k27me3.txt"))
  dat.meta <- dat.meta.lst[[jmark]] %>%
    dplyr::select(cell, ctype.from.LL, batch) %>%
    rowwise() %>%
    mutate(cluster = ctype.from.LL,
           jrep2 = ifelse(batch == "New", "anew", "old"),
           cluster = ifelse(cluster == "HSCs", "aHSCs", cluster))
  
  dat.meta$ctype.from.LL <- NULL
  dat.meta$batch <- NULL
  
  
  
  # Go genome-wide ----------------------------------------------------------
  
  rnames.keep <- rownames(dat.impute.log2)  # keep all
  jmat.long <- dat.impute.log2[rnames.keep, ] %>%
    data.table::melt()
  colnames(jmat.long) <- c("rname", "cell", "log2exprs")
  jmat.long <- jmat.long %>%
    left_join(., dat.meta)
  
  print("Correcting batch effects... multicore")
  
  # multicore
  jmat.long.lst <- split(jmat.long, jmat.long$rname)
  
  # clean up 
  rm(jmat.long, dat.impute.log2, out.lda)
  
  jmat.wide.adj.lst <- parallel::mclapply(jmat.long.lst, function(jmat.long){
    jmat.long.adj <- jmat.long %>%
      group_by(rname) %>%
      do(AdjustBatchEffect2(.))
    
    # make wide vector
    jmat <- jmat.long.adj %>%
      reshape2::dcast(data = ., formula = "rname ~ cell", value.var = "log2exprsadj")
    return(jmat)
    # create 
  }, mc.cores = 8) %>%
    bind_rows()
  
  # jmat.wide.adj.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_long_k27me3_batch_corrected.2022-04-19.rds")
  saveRDS(jmat.wide.adj.lst, file = outrds.lst)
  
  jmat.wide <- jmat.wide.adj.lst %>%
    bind_rows()
  rnames.keep <- jmat.wide$rname
  jmat.wide$rname <- NULL
  jmat.wide <- as.matrix(jmat.wide)
  rownames(jmat.wide) <- rnames.keep
  
  # arrange colnames to match
  cnames <- colnames(count.mat)
  rnames <- rownames(count.mat)
  jmat.wide <- jmat.wide[rnames, cnames]
  
  assertthat::assert_that(all(rownames(count.mat) == rownames(jmat.wide)))
  assertthat::assert_that(all(colnames(count.mat) == colnames(jmat.wide)))
  
  saveRDS(jmat.wide, file = outrds)
  
  # save raw counts
  saveRDS(count.mat, file = outrdsraw)
  
  print("Correcting batch effects... done")
  print(Sys.time() - jstart)
  
  
  
}



