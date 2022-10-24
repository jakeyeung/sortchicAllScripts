# Jake Yeung
# Date of Creation: 2022-05-03
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/28-pseudotime_cubic_spline_each_gene_get_derivatives.R
# 

rm(list=ls())

library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(mgcv)
library(gratia)

jratio <- 0.66

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks



# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta)
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

# Load LDAs: all  ---------------------------------------------------------


tm.lst <- lapply(jmarks, function(jmark){
  if (jmark == "k4me1"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    # inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-15/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-15.Robj")
    # inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_new_only_dynbins.k27me3.2022-04-15/ldaOut.count_mat_new_only.", jmark, ".2022-04-15.Robj")
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_new_only_dynbins.k27me3.2022-04-15/ldaOut.count_mat_new_only_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  
  load(inf.ldaout, v=T)
  tm <- posterior(out.lda)
  return(tm)
})

# assertthat::assert_that(nrow(dat.meta.lst$k27me3) == nrow(tm.lst$k27me3$topics))

print(lapply(tm.lst, function(x) dim(x$topics)))

dat.impute.lst <- lapply(tm.lst, function(tm){
  dat.impute <- log2(t(tm$topics %*% tm$terms))
})


# Load batch corrected k27me3 ---------------------------------------------------

inf.impute.k27me3.bc <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_batch_corrected.2022-04-19.rds"
dat.impute.k27me3.bc <- readRDS(inf.impute.k27me3.bc)
dat.impute.lst$k27me3 <- dat.impute.k27me3.bc

# Load trajs --------------------------------------------------------------

# indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix"
indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
dat.trajs <- lapply(jmarks, function(jmark){
  inf.trajs <- file.path(indir.traj, paste0("trajs_outputs.", jmark, ".rds"))
  if (jmark == "k27me3"){
    inf.trajs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix/trajs_outputs_batch_corrected.k27me3.2022-04-21.rds")
  }
  print(inf.trajs)
  readRDS(inf.trajs)
})

dat.trajs$k27me3$Granulocytes <- subset(dat.trajs$k27me3$Granulocytes, ctype.from.LL != "Basophils")


# Fit gam for each region  --------------------------------------------------

# jmark <- "k27me3"
# set.seed(0)
# # dat.impute.lst$k27me3
# jregion <- sample(rownames(dat.impute.lst[[jmark]]), size = 1)
# print(jregion)
outrds.final <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/gam_fits_dynamic_bins4.final.", Sys.Date(), ".rds")


# fits.lst.bymark <- parallel::mclapply(jmarks, function(jmark){
fits.lst.bymark <- lapply(jmarks, function(jmark){
  outrds.tmp <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/gam_fits_dynamic_bins4.", jmark, ".", Sys.Date(), ".rds")
  jregions <- rownames(dat.impute.lst[[jmark]]); names(jregions) <- jregions
  
  jfits.lst.byregion <- lapply(jregions, function(jregion){
    
    dat.signal <- data.frame(cell = colnames(dat.impute.lst[[jmark]]), signal = dat.impute.lst[[jmark]][jregion, ], stringsAsFactors = FALSE)
    
    traj.ctypes <- names(dat.trajs[[jmark]])
    names(traj.ctypes) <- traj.ctypes
    dat.trajs.long <- lapply(traj.ctypes, function(traj.ctype){
      dat.trajs.sub <- dat.trajs[[jmark]][[traj.ctype]] %>%
        filter(is.ctype) %>%
        dplyr::select(cell, ctype.ordered, ptime) %>%
        mutate(traj = traj.ctype)
      # colcode = ctype2col[[traj.ctype]])
    }) %>%
      bind_rows() %>%
      left_join(., dat.signal, by = "cell") %>%
      left_join(., subset(dat.meta.lst[[jmark]], select = c(cell, batch, ctype.from.LL, colcode)), by = "cell") %>%
      mutate(traj = as.factor(traj),
             region = jregion)
    
    dat.trajs.long.split <- split(dat.trajs.long, dat.trajs.long$traj)
    
    jfits.split <- lapply(dat.trajs.long.split, function(jsub){
      jfit.sub <- gam(formula = signal ~ s(ptime, k = 4, bs = "cs", by = traj), gamma = 10, method = "REML", data = jsub)
      jsub$pred <- predict(jfit.sub)
      # derivatives(jfit.sub)
      return(list(jsub = jsub, fit = jfit.sub))
    })
    return(jfits.split)
  })
  print(paste("Done for", jmark))
  saveRDS(jfits.lst.byregion, file = outrds.tmp)
  return(jfits.lst.byregion)
# }, mc.cores = length(jmarks))
})

print(fits.lst.bymark)
saveRDS(fits.lst.bymark, file = outrds.final)
