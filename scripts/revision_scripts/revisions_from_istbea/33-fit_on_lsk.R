# Jake Yeung
# Date of Creation: 2022-05-06
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/33-fit_on_lsk.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# Plot LSK  ---------------------------------------------------------------

dat.sub.impute.knn.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat.sub.impute.knn.lst.2022-05-06.rds")

jmarks <- names(dat.sub.impute.knn.lst); names(jmarks) <- jmarks

dat.lsk.long <- lapply(jmarks, function(jmark){
  jsub <- dat.sub.impute.knn.lst[[jmark]]
  jsub$mark <- jmark
  return(jsub)
}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(ptime.lsk.unnorm = -1 * sca1_f + 1 * lin_f) %>%
  ungroup() %>%
  mutate(ptime.lsk = ( ptime.lsk.unnorm - (min(ptime.lsk.unnorm)) ) / (max(ptime.lsk.unnorm) - min(ptime.lsk.unnorm)) )

outdir.rds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs"
cell2ptime.lsk <- hash::hash(dat.lsk.long$cell, dat.lsk.long$ptime.lsk)
saveRDS(dat.lsk.long, file.path(outdir.rds, paste0("dat_lsk_long_ptime_lsk.", Sys.Date(), ".rds")))
saveRDS(cell2ptime.lsk, file.path(outdir.rds, paste0("cell2ptime_lsk.", Sys.Date(), ".rds")))


ggplot(dat.lsk.long, aes(x = ptime.lsk)) + 
  geom_histogram() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

lsk.mat <- as.matrix(subset(dat.lsk.long, select = c(sca1_f, ckit_f, lin_f)))
rownames(lsk.mat) <- dat.lsk.long$cell


# Load exprs mats ---------------------------------------------------------------




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
  dat.trajs.bytraj <- readRDS(inf.trajs)
  dat.trajs.bytraj <- lapply(dat.trajs.bytraj, function(jdat){
    jdat %>%
      rowwise() %>%
      mutate(ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
      filter(!is.na(ptime.lsk))
  })
})

dat.trajs$k27me3$Granulocytes <- subset(dat.trajs$k27me3$Granulocytes, ctype.from.LL != "Basophils")


# Fit gam for each region  --------------------------------------------------

outrds.final <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/gam_fits_dynamic_bins_on_LSK_bins.", Sys.Date(), ".rds")

set.seed(0)
fits.lst.bymark <- lapply(jmarks, function(jmark){
  outrds.tmp <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/gam_fits_dynamic_bins_on_LSK_bins.", jmark, ".", Sys.Date(), ".rds")
  jregions <- rownames(dat.impute.lst[[jmark]]); names(jregions) <- jregions
  
  # (jregion <- sample(jregions, size = 1))
  # jregions <- sample(x = jregions, size = round(length(jregions) / 5), replace = FALSE)
  jfits.lst.byregion <- lapply(jregions, function(jregion){
    
    dat.signal <- data.frame(cell = colnames(dat.impute.lst[[jmark]]), signal = dat.impute.lst[[jmark]][jregion, ], stringsAsFactors = FALSE)
    
    traj.ctypes <- names(dat.trajs[[jmark]])
    names(traj.ctypes) <- traj.ctypes
    dat.trajs.long <- lapply(traj.ctypes, function(traj.ctype){
      dat.trajs.sub <- dat.trajs[[jmark]][[traj.ctype]] %>%
        filter(is.ctype) %>%
        dplyr::select(cell, ctype.ordered, ptime.lsk) %>%
        mutate(traj = traj.ctype)
      # colcode = ctype2col[[traj.ctype]])
    }) %>%
      bind_rows() %>%
      left_join(., dat.signal, by = "cell") %>%
      # left_join(., subset(dat.meta.lst[[jmark]], select = c(cell, batch, ctype.from.LL, colcode)), by = "cell") %>%
      mutate(region = jregion) %>%
      mutate(ptime.lsk.bin = floor(ptime.lsk / 0.05) * 0.05)
    
    # bin by traj
    
    dat.trajs.long.split <- split(dat.trajs.long, dat.trajs.long$traj)
    
    jfits.split <- lapply(dat.trajs.long.split, function(jsub){
      jsub.traj <- jsub %>%
        group_by(ptime.lsk.bin, traj) %>%
        summarise(signal.withintraj.mad = mean(abs(signal - mean(signal))),
                  signal.withintraj.sd = mean( (signal - mean(signal) ) ^ 2 ), 
                  signal.mean = mean(signal),
                  .groups = "drop")
      jsub.global <- subset(jsub, ptime.lsk < 0.25) %>%
        ungroup() %>%
        summarise(signal.global.mad = mean(abs(signal - mean(signal))),
                  signal.global.sd = mean( (signal - mean(signal) ) ^ 2 ), 
                  signal.global.mean = mean(signal),
                  .groups = "drop")
      
      return(list(jsub.traj = jsub.traj, jsub.global = jsub.global))
        
      # jfit.sub <- gam(formula = signal ~ s(ptime.lsk, k = 4, bs = "cs", by = traj), gamma = 10, method = "REML", data = jsub)
      # jsub$pred <- predict(jfit.sub)
      # # derivatives(jfit.sub)
      # return(list(jsub = jsub, fit = jfit.sub))
    })
    
    
    # # plot outputs
    # jcheck <- lapply(jfits.split, function(jsub){
    #   jsub$jsub.traj
    # }) %>%
    #   bind_rows()
    # 
    # jcheck.sum <- jcheck %>%
    #   group_by(ptime.lsk.bin) %>%
    #   summarise(signal.mad = mean ( abs(signal.mean - mean(signal.mean) )), 
    #             signal.sd = mean ( abs ( signal.mean - mean(signal.mean) ) ^ 2 ))
    
    # ggplot() + 
    #   geom_point(mapping = aes(x = ptime.lsk, y = signal, color = traj), data = dat.trajs.long, alpha = 0.05) + 
    #   geom_line(mapping = aes(x = ptime.lsk.bin, y = signal.mean, color = traj, group = traj), data = jcheck) + 
    #   theme_bw() + 
    #   ggtitle(jregion) + 
    #   # facet_wrap(~traj) + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    # ggplot() + 
    #   geom_line(mapping = aes(x = ptime.lsk.bin, y = signal.mad), data = jcheck.sum) + 
    #   theme_bw() + 
    #   ggtitle(jregion) + 
    #   # facet_wrap(~traj) + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    # # ggplot(jcheck, aes(x = ptime.lsk.bin, y = pred, color = traj)) + 
    # #   geom_point() + 
    # #   theme_bw() + 
    # #   ggtitle(jregion) + 
    # #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    return(jfits.split)
  })
  print(paste("Done for", jmark))
  
  # check downstream
  jfits.lst.long <- lapply(jfits.lst.byregion, function(jfit.byctype){
    lapply(jfit.byctype, function(jfit){
      jfit$jsub.traj
    }) %>%
      bind_rows() %>%
      group_by(ptime.lsk.bin) %>%
      summarise(signal.mad = mean ( abs(signal.mean - mean(signal.mean) )), 
                signal.sd = mean ( abs ( signal.mean - mean(signal.mean) ) ^ 2 ))
      
  }) %>%
    bind_rows()
  
  jfits.global.long <- lapply(jfits.lst.byregion, function(jfit.byctype){
    lapply(jfit.byctype, function(jfit){
      jfit$jsub.global
    }) %>%
      bind_rows() 
  }) %>%
    bind_rows() %>%
    ungroup() %>%
    summarise(signal.global.mean = mean(signal.global.mean), 
              signal.global.mad = mean(signal.global.mad), 
              signal.global.sd = mean(signal.global.sd))
  
  jfits.lst.long$signal.global.mean <- jfits.global.long$signal.global.mean
  jfits.lst.long$signal.global.mad <- jfits.global.long$signal.global.mad
  jfits.lst.long$signal.global.sd <- jfits.global.long$signal.global.sd
  
  ggplot() + 
    geom_boxplot(mapping = aes(x = ptime.lsk.bin, y = signal.mad, group = ptime.lsk.bin), data = jfits.lst.long) + 
    geom_hline(yintercept = 0.3) + 
    # geom_ribbon(mapping = aes(x = ptime.lsk.bin, ymin = 0, ymax = signal.global.mad, group = ptime.lsk.bin), data = jfits.lst.long) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
  
  saveRDS(jfits.lst.byregion, file = outrds.tmp)
  return(jfits.lst.byregion)
})




# For each gene, for each mark, for each trajectory, fit GAM across pseudotime 




# Calculate deviation from mean at every binned pseudotime  ---------------












