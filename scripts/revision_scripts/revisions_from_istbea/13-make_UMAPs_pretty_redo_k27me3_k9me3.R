# Jake Yeung
# Date of Creation: 2022-02-15
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/13-make_UMAPs_pretty_redo_k27me3_k9me3.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks


# Load metas --------------------------------------------------------------

indir.meta.init <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots"
dat.meta.lst1 <- lapply(jmarks, function(jmark){
  fread(file.path(indir.meta.init, paste0("metadata_with_colors.", jmark, ".2022-02-13.txt")))
})


# Load  metas for repressive -------------------------------------------





# Replace metas for some marks --------------------------------------------





# Plot outputs ------------------------------------------------------------





# Show batch effects ------------------------------------------------------






# plot pretty metas  ------------------------------------------------------


