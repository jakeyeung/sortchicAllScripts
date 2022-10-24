# Jake Yeung
# Date of Creation: 2022-07-26
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/8-make_metadata_for_splitting_bams_by_celltype.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks

inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.clean.lst <- lapply(jmarks, function(jmark){
  dat.meta.lst[[jmark]] %>%
    dplyr::select(cell, ctype.from.LL)
})



# Write cells in column 1, celltype in column 2  --------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata"

new2old <- hash::hash(jmarks, jmarksold)

for (jmark in jmarks){
  print(jmark)
  jmarkold <- new2old[[jmark]]
  fwrite(dat.meta.clean.lst[[jmark]], file = file.path(outdir, paste0("cell_ctype_metadata_for_splitting_bams.", jmarkold, ".txt")), sep = "\t", col.names = FALSE)
}

