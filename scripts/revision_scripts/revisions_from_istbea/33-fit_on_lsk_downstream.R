# Jake Yeung
# Date of Creation: 2022-05-07
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/33-fit_on_lsk_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

jfits.lst.bymark <- lapply(jmarks, function(jmark){
  inf.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_lsk/gam_fits_lsk_dynamic_bins.", jmark, ".fromargs.full_list.summarize_mean_mad_sd.2022-05-07.rds")
  readRDS(inf.out)
})

jfits.lst.bymark.summary <- lapply(jfits.lst.bymark, function(jfits.lst.byregion){
  
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
  
  return(jfits.lst.long)
})


m.lst <- lapply(jmarks, function(jmark){
  jfits <- jfits.lst.bymark.summary[[jmark]]
  ggplot() + 
    geom_boxplot(mapping = aes(x = ptime.lsk.bin, y = signal.mad, group = ptime.lsk.bin), data = jfits) + 
    geom_hline(mapping = aes(yintercept = signal.global.mad), data = jfits) + 
    # geom_ribbon(mapping = aes(x = ptime.lsk.bin, ymin = 0, ymax = signal.global.mad, group = ptime.lsk.bin), data = jfits.lst.long) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst)
