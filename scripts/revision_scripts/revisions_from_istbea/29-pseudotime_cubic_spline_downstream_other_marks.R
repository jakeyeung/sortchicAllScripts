# Jake Yeung
# Date of Creation: 2022-05-04
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/29-pseudotime_cubic_spline_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(mgcv)
library(gratia)
library(scchicFuncs)
library(topicmodels)


deltaT <- 0.01

jmarksall <- c("k4me3", "k9me3", "k27me3", "k4me1"); names(jmarksall) <- jmarksall
jmarks <- c("k27me3", "k9me3", "k4me1"); names(jmarks) <- jmarks
# jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks



# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarksall, function(jmark){
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



# Load cubic spline  ------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  
  outf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/signal_extrap_vectors_regions_ctypes_deltaT.", deltaT, ".", jmark, ".", Sys.Date(), ".rds")
  # dat.fits <- readRDS(paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/gam_fits_dynamic_bins.", jmark, ".2022-05-04.rds"))
  dat.fits <- readRDS(paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/gam_fits_dynamic_bins.", jmark, ".fromargs.full_list.rds"))
  
  # Get derivatives ---------------------------------------------------------
  
  # plot fit
  
  ctypes.all <- names(dat.trajs[[jmark]]); names(ctypes.all) <- ctypes.all
  jregions <- names(dat.fits); names(jregions) <- jregions
  
  
  # jregion <- "chr1:9950000-10000000"
  # jctype <- "Eryths"
  
  
  jdat.merge.lst <- lapply(jregions, function(jregion){
    jdat.merge.byctype <- lapply(ctypes.all, function(jctype){
      
      
      jdat <- dat.fits[[jregion]][[jctype]]$jsub %>%
        arrange(ptime)
      jdat$pred <- predict(dat.fits[[jregion]][[jctype]]$fit)
      # print(head(jdat))
      jderiv <- derivatives(dat.fits[[jregion]][[jctype]]$fit, newdata = jdat)
      
      jderiv.clean <- jderiv %>%
        dplyr::select(data, derivative) %>%
        dplyr::rename(ptime = data)
      
      jdat.merge <- left_join(jdat, jderiv.clean, by = "ptime")
      
      jdat.merge <- jdat.merge %>%
        mutate(signal.extrap = signal + derivative * deltaT)
      return(jdat.merge)
    })  %>%
      bind_rows()
    jdat.merge.byctype.sum <- jdat.merge.byctype %>%
      group_by(cell) %>%
      summarise(signal.extrap = mean(signal.extrap), 
                derivative = mean(derivative), 
                pred = mean(pred), 
                ptime = mean(ptime), 
                signal = unique(signal))
    # return signal and extrap vector
    signal.vec <- jdat.merge.byctype.sum$signal
    signal.extrap.vec <- jdat.merge.byctype.sum$signal.extrap
    
    names(signal.vec) <- jdat.merge.byctype.sum$cell
    names(signal.extrap.vec) <- jdat.merge.byctype.sum$cell
    
    return(list(signal.vec = signal.vec, signal.extrap.vec = signal.extrap.vec))
  })
  saveRDS(jdat.merge.lst, file = outf)
}



    
#     
# m1 <- ggplot() +
#   geom_point(mapping = aes(ptime, y = signal), data = jdat, alpha = 0.1, color = "grey85") + 
#   geom_path(mapping = aes(ptime, y = pred), data = jdat, color = "blue") + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# m2 <- ggplot(jderiv, aes(x = data, y = derivative)) + 
#   geom_line() + 
#   geom_hline(yintercept = 0, linetype = "dotted") + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# JFuncs::multiplot(m1, m2, cols = 1)
# 
# jdat.merge <- left_join(jdat, jderiv.clean, by = "ptime")
# 
# 
# # Predict next timepoints from deriv --------------------------------------
# 
# deltaT <- 0.01
# 
# jdat.merge <- jdat.merge %>%
#   mutate(signal.extrap = signal + derivative * deltaT)
# 
# 
# # head(derivatives(dat.fits$`chr1:9950000-10000000`$Eryths$fit))
# 
# 
# # Plot PCA with extrap ----------------------------------------------------


