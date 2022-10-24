# Jake Yeung
# Date of Creation: 2022-08-02
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/12-filter_out_SLC5_CMP.R
# 
rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks


# load --------------------------------------------------------------------

indir.mat <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_progenitors_filt"

count.mat.tss.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(indir.mat, paste0("countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27.rds"))
  print(inf)
  count.mat <- readRDS(inf)
  return(count.mat)
})

count.mat.dynbins.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(indir.mat, paste0("countmat_dynbins_new_progenitors_only.", jmark, ".2022-07-27.rds"))
  count.mat <- readRDS(inf)
  return(count.mat)
})


cells.prog.lst <- lapply(count.mat.tss.lst, function(jmat){
  colnames(jmat)
})


cells.prog.nosl5.lst <- lapply(cells.prog.lst, function(cells.vec){
  cells.vec[!grepl("SL5", cells.vec)]
})

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

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

cells.prog.nosl5.nocmp.lst <- lapply(jmarks, function(jmark){
  cells.cmp <- subset(dat.meta.lst[[jmark]], ctype.from.LL == "CMP")$cell
  cells.vec <- cells.prog.nosl5.lst[[jmark]]
  cells.filt.vec <- cells.vec[!cells.vec %in% cells.cmp]
  return(cells.filt.vec)
})

# Filter out SLC5 ---------------------------------------------------------

count.mat.TSS.nosl5.lst <- lapply(jmarks, function(jmark){
  mat <- count.mat.tss.lst[[jmark]]
  cells.keep <- cells.prog.nosl5.lst[[jmark]]
  mat.filt <- mat[, cells.keep]
})


count.mat.dynbins.nosl5.lst <- lapply(jmarks, function(jmark){
  mat <- count.mat.dynbins.lst[[jmark]]
  cells.keep <- cells.prog.nosl5.lst[[jmark]]
  mat.filt <- mat[, cells.keep]
})



# Filter out CMP additionally ---------------------------------------------

count.mat.TSS.nosl5.nocmp.lst <- lapply(jmarks, function(jmark){
  mat <- count.mat.tss.lst[[jmark]]
  cells.keep <- cells.prog.nosl5.nocmp.lst[[jmark]]
  mat.filt <- mat[, cells.keep]
})

count.mat.dynbins.nosl5.nocmp.lst <- lapply(jmarks, function(jmark){
  mat <- count.mat.dynbins.lst[[jmark]]
  cells.keep <- cells.prog.nosl5.nocmp.lst[[jmark]]
  mat.filt <- mat[, cells.keep]
})

# Write outputs -----------------------------------------------------------


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_progenitors_filt_further"
dir.create(outdir)
# outdir.nosl5.nocmp <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_progenitors_filt_nosl5_nocmp"

for (jmark in jmarks){
  print(jmark)
  outf.tss.nosl5 <- file.path(outdir, paste0("countmat_TSS_new_progenitors_only.nosl5.", jmark, ".", Sys.Date(), ".rds"))
  outf.dynbins.nosl5 <- file.path(outdir, paste0("countmat_dynbins_new_progenitors_only.nosl5.", jmark, ".", Sys.Date(), ".rds"))
  
  outf.tss.nosl5.nocmp <- file.path(outdir, paste0("countmat_TSS_new_progenitors_only.nosl5.nocmp", jmark, ".", Sys.Date(), ".rds"))
  outf.dynbins.nosl5.nocmp <- file.path(outdir, paste0("countmat_dynbins_new_progenitors_only.nosl5.nocmp", jmark, ".", Sys.Date(), ".rds"))
  
  count.mat.TSS.nosl5.tmp <- count.mat.TSS.nosl5.lst[[jmark]]
  count.mat.TSS.nosl5.nocmp.tmp <- count.mat.TSS.nosl5.nocmp.lst[[jmark]]
  
  count.mat.dynbins.nosl5.tmp <- count.mat.dynbins.nosl5.lst[[jmark]]
  count.mat.dynbins.nosl5.nocmp.tmp <- count.mat.dynbins.nosl5.nocmp.lst[[jmark]]
  
  print(dim(count.mat.TSS.nosl5.tmp))
  print(dim(count.mat.TSS.nosl5.nocmp.tmp))
  print(dim(count.mat.dynbins.nosl5.tmp))
  print(dim(count.mat.dynbins.nosl5.nocmp.tmp))
  
  
  saveRDS(count.mat.TSS.nosl5.tmp, file = outf.tss.nosl5)
  saveRDS(count.mat.dynbins.nosl5.tmp, file = outf.dynbins.nosl5)
  
  saveRDS(count.mat.TSS.nosl5.nocmp.tmp, file = outf.tss.nosl5.nocmp)
  saveRDS(count.mat.dynbins.nosl5.nocmp.tmp, file = outf.dynbins.nosl5.nocmp)
}

