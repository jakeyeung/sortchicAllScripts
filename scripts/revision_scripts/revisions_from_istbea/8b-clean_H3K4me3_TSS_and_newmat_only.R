# Jake Yeung
# Date of Creation: 2022-04-12
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8b-clean_H3K4me3_TSS_and_newmat_only.R
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

inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS/ldaAnalysis_fripfilt_TSS/ldaOut.count_mat_TSS_combined.", jmark, ".2022-02-07.Robj")
load(inf.lda, v=T)
count.mat.tss <- count.mat

indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
inf.ldaout <- file.path(indir.ldaout, paste0("ldaOut.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj"))
load(inf.ldaout, v=T)

count.mat.dynbins <- count.mat


# Load clean meta ---------------------------------------------------------

inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata_k4me3_cleaned_up/metadata_k4me3_louvain_filtered.no3no7.2022-04-12.txt"
dat.meta <- fread(inf.meta)

cells.keep <- dat.meta$cell
length(cells.keep)
assertthat::assert_that(unique(length(cells.keep) == length(cells.keep)))

cells.keep.new <- subset(dat.meta, batch == "New")$cell
length(cells.keep.new)
assertthat::assert_that(unique(length(cells.keep.new) == length(cells.keep.new)))

# Save output -------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k4me3_remove_bad_louvains_dynbins_TSS"
dir.create(outdir)

fname.dynbins.all <- paste0("count_mat_cleaned_no3no7_dynbins_allcells.",  jmark, ".", Sys.Date(), ".rds")
fname.dynbins.new <- paste0("count_mat_cleaned_no3no7_dynbins_newcells.",  jmark, ".", Sys.Date(), ".rds")
fname.tss.all <- paste0("count_mat_cleaned_no3no7_tss_allcells.",  jmark, ".", Sys.Date(), ".rds")
fname.tss.new <- paste0("count_mat_cleaned_no3no7_tss_newcells.",  jmark, ".", Sys.Date(), ".rds")

saveRDS(count.mat.dynbins[, cells.keep], file = file.path(outdir, fname.dynbins.all))
saveRDS(count.mat.dynbins[, cells.keep.new], file = file.path(outdir, fname.dynbins.new))
saveRDS(count.mat.tss[, cells.keep], file = file.path(outdir, fname.tss.all))
saveRDS(count.mat.tss[, cells.keep.new], file = file.path(outdir, fname.tss.new))




