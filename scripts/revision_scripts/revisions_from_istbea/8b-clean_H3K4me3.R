# Jake Yeung
# Date of Creation: 2022-04-11
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8b-clean_H3K4me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Clean up k4me3  ---------------------------------------------------------

jsuffix <- "dynamicbins"

jmark <- "k4me3"
indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
inf.ldaout <- file.path(indir.ldaout, paste0("ldaOut.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj"))
load(inf.ldaout, v=T)


indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap"
inf.meta <- file.path(indir.meta, paste0("metadata.", jmark, ".2022-03-01.txt"))
dat.meta <- fread(inf.meta) %>%
  dplyr::select(-ctype.from.LL)


# Remove louvains ---------------------------------------------------------

bad.louvains <- c("3", "7")

dat.meta.no3 <- subset(dat.meta, louvain != "3")
dat.meta.no7 <- subset(dat.meta, louvain != "7")
dat.meta.no3no7 <- subset(dat.meta, !louvain %in% c("3", "7"))

print(dim(dat.meta.no3))
print(dim(dat.meta.no7))
print(dim(dat.meta.no3no7))

cells.keep.no3 <- subset(dat.meta.no3)$cell
cells.keep.no7 <- subset(dat.meta.no7)$cell
cells.keep.no3no7 <- subset(dat.meta.no3no7)$cell

count.mat.no3 <- count.mat[, cells.keep.no3]
count.mat.no7 <- count.mat[, cells.keep.no7]
count.mat.no3no7 <- count.mat[, cells.keep.no3no7]

print(dim(count.mat.no3))
print(dim(count.mat.no7))
print(dim(count.mat.no3no7))

# Save output -------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k4me3_remove_bad_louvains"
fname.no3 <- paste0("count_mat_cleaned_no3.", jmark, ".", Sys.Date(), ".rds")
fname.no7 <- paste0("count_mat_cleaned_no7.", jmark, ".", Sys.Date(), ".rds")
fname.no3no7 <- paste0("count_mat_cleaned_no3no7.", jmark, ".", Sys.Date(), ".rds")

saveRDS(count.mat.no3, file = file.path(outdir, fname.no3))
saveRDS(count.mat.no7, file = file.path(outdir, fname.no7))
saveRDS(count.mat.no3no7, file = file.path(outdir, fname.no3no7))


