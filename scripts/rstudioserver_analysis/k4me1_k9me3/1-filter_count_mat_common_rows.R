# Jake Yeung
# Date of Creation: 2020-10-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/1-filter_count_mat_common_rows.R
# Set up for LDA 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"

# 
# # Load singles from csv ---------------------------------------------------
# 
# indir.raw <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/countTablesAndRZr1only_ByChromo.NewFilters.bybed.blfix")
# assertthat::assert_that(dir.exists(indir.raw))
# infs.k4me1 <- list.files(indir.raw, pattern = ".*.K4me1.*dbl_bins.csv", full.names = TRUE)
# infs.k9me3 <- list.files(indir.raw, pattern = ".*.K9me3.*dbl_bins.csv", full.names = TRUE)
# 
# mat.k4me1.raw.lst <- lapply(infs.k4me1, function(inf) ReadMatTSSFormat(inf))
# k4me1.common.rows <- Reduce(f = intersect, x = lapply(mat.k4me1.raw.lst, function(x) rownames(x)))
# 
# mat.k4me1.raw.merged.lst <- lapply(mat.k4me1.raw.lst, function(jmat){
#   jmat[k4me1.common.rows, ]
# })
# 
# 
# mat.k9me3.raw.lst <- lapply(infs.k9me3, function(inf) ReadMatTSSFormat(inf))
# k9me3.common.rows <- Reduce(f = intersect, x = mat.k9me3.raw.lst)
# 
# mat.k4me1.raw <- do.call(cbind, mat.k4me1.raw.lst)
# mat.k9me3.raw <- do.call(cbind, mat.k9me3.raw.lst)

# Load singles from filtered ------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/varfilt/count_mat.H3K4me1.varfilt_1.5.rds"))
inf.k9me3 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/varfilt/count_mat.H3K9me3.varfilt_1.5.rds"))

mat.k4me1 <- readRDS(inf.k4me1)
mat.k9me3 <- readRDS(inf.k9me3)

good.cells.k4me1 <- colnames(mat.k4me1)
good.cells.k9me3 <- colnames(mat.k9me3)

bins.k4me1 <- rownames(mat.k4me1)
bins.k9me3 <- rownames(mat.k9me3)

# Get good bins -----------------------------------------------------------

inf.dbl <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl.blfix/count_mat.H3K9me3-H3K4me1.filt_0.15_0.95_counts_and_l2r.blfix.rds")
mat.dbl <- readRDS(inf.dbl)

cols.keep.dbl <- colnames(mat.dbl)

bins.dbl <- rownames(mat.dbl)

bins.lst <-  list(bins.k4me1, bins.k9me3, bins.dbl)

bins.keep <- Reduce(f = intersect, x = bins.lst)
# bins.keep <- dat.bed$V4

lapply(bins.lst, length)

length(bins.keep)


# Get common bins across all ----------------------------------------------




# Write filterd mats ------------------------------------------------------

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.single_match_dbl")

outf1 <- file.path(outdir, "count_mat.H3K4me1.match_dbl.rds")
outf2 <- file.path(outdir, "count_mat.H3K9me3.match_dbl.rds")
outf3 <- file.path(outdir, "count_mat.H3K4me1xH3K9me3.match_dbl.rds")


cols.keep.k4me1 <- colnames(mat.k4me1) %in% good.cells.k4me1
cols.keep.k9me3 <- colnames(mat.k9me3) %in% good.cells.k9me3

rows.keep.k4me1 <- rownames(mat.k4me1) %in% bins.keep
rows.keep.k9me3 <- rownames(mat.k9me3) %in% bins.keep

mat.k4me1.filt <- mat.k4me1[bins.keep, cols.keep.k4me1]
mat.k9me3.filt <- mat.k9me3[bins.keep, cols.keep.k9me3]
mat.dbl.filt <- mat.dbl[bins.keep, cols.keep.dbl]

saveRDS(mat.k4me1.filt, file = outf1)
saveRDS(mat.k9me3.filt, file = outf2)
saveRDS(mat.dbl.filt, file = outf3)




