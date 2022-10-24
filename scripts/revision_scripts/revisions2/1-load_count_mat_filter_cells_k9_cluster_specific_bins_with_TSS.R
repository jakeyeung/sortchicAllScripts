# Jake Yeung
# Date of Creation: 2022-07-27
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/1-load_count_mat_filter_cells_k9_cluster_specific_bins_with_TSS.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

keeptop <- 500

metadir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/processed_data"
# matmain <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_k9_cluster_specific_bins_keeptop_", keeptop, "/counts_tables_k9_cluster_specific_bins")
matmain <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_k9_cluster_specific_bins_keeptop_", keeptop, "_with_TSS/counts_tables_k9_cluster_specific_bins")
assertthat::assert_that(dir.exists(matmain))
outdir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/LDA_count_mats_inputs_k9_cluster_specific_bins_keeptop_", keeptop, "_with_TSS")
dir.create(outdir)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew


# Load metadata -----------------------------------------------------------

# dat.metas <- lapply(jmarksnew, function(jmark){
#   print(jmark)
#   fread(file.path(metadir, paste0("metadata_bonemarrow_allmerged_", jmark, ".txt.gz")))
# })


indir.meta.bm.clean <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed"

# dat.metas.lsk.bm <- lapply(jmarksnew, function(jmark){
#   inf.meta.bm <- file.path(indir.meta.bm.clean, paste0("LSK_metadata_colorcode_more_contrast.", jmark, ".2022-05-17.txt"))
#   fread(inf.meta.bm)
# })

dat.metas <- lapply(jmarksnew, function(jmark){
  inf.meta.bm <- file.path(indir.meta.bm.clean, paste0("umap_metadata_color_DC_monocyte_fixed.", jmark, ".2022-05-17.txt"))
  fread(inf.meta.bm)
})




# Load mats ---------------------------------------------------------------

dat.mats <- lapply(jmarks, function(jmark){
  print(jmark)
  matdir <- file.path(matmain, paste0("BM_", jmark))
  inf <- file.path(matdir, paste0("BM_allmerged_", jmark, ".countTable.binsize_k9_cluster_specific_bins.csv"))
  print(inf)
  ReadMatTSSFormat(inf = inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = TRUE)
})


# Filter mats -------------------------------------------------------------

names(dat.metas) <- jmarks
dat.mats.filt <- lapply(jmarks, function(jmark){
  print(jmark)
  cells.keep <- dat.metas[[jmark]]$cell
  dat.mat.tmp <- dat.mats[[jmark]]
  cells.keep2 <- colnames(dat.mat.tmp) %in% cells.keep
  dat.mat.filt <- dat.mat.tmp[, cells.keep2]
})


rnames.cleaned.lst <- lapply(dat.mats.filt, function(jdat) rownames(jdat))
# rnames.cleaned.lst <- rnames.lst
rnames.cleaned.common <- Reduce(f = intersect, x = rnames.cleaned.lst)
print(sapply(rnames.cleaned.lst, length))
print(length(rnames.cleaned.common))

dat.mats.filt.rnames.cleaned <- lapply(jmarks, function(jmark){
  mat.tmp <- dat.mats.filt[[jmark]]
  rownames(mat.tmp) <- rnames.cleaned.lst[[jmark]]
  # get common
  mat.tmp <- mat.tmp[rnames.cleaned.common, ]
  return(mat.tmp)
})

lapply(dat.mats.filt.rnames.cleaned, dim)



# Write output ------------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outftmp <- file.path(outdir, paste0("count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_", keeptop, "_with_TSS.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(object = dat.mats.filt.rnames.cleaned[[jmark]], file = outftmp)
}




