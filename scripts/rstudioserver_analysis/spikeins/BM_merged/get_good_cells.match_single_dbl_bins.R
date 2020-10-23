# Jake Yeung
# Date of Creation: 2020-10-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/get_good_cells.R
# Filter good cells

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load singles from csv ---------------------------------------------------

indir.raw <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/countTablesAndRZr1only_ByChromo.NewFilters.bybed.blfix")
assertthat::assert_that(dir.exists(indir.raw))
infs.k4me1 <- list.files(indir.raw, pattern = ".*.K4me1.*dbl_bins.csv", full.names = TRUE)
infs.k9me3 <- list.files(indir.raw, pattern = ".*.K9me3.*dbl_bins.csv", full.names = TRUE)

mat.k4me1.raw.lst <- lapply(infs.k4me1, function(inf) ReadMatTSSFormat(inf))
k4me1.common.rows <- Reduce(f = intersect, x = lapply(mat.k4me1.raw.lst, function(x) rownames(x)))

mat.k4me1.raw.merged.lst <- lapply(mat.k4me1.raw.lst, function(jmat){
  jmat[k4me1.common.rows, ]
})


mat.k9me3.raw.lst <- lapply(infs.k9me3, function(inf) ReadMatTSSFormat(inf))
k9me3.common.rows <- Reduce(f = intersect, x = mat.k9me3.raw.lst)

mat.k4me1.raw <- do.call(cbind, mat.k4me1.raw.lst)
mat.k9me3.raw <- do.call(cbind, mat.k9me3.raw.lst)

# Load singles from filtered ------------------------------------------------------------

inf.k4me1 <- file.path(hubprefix, paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/varfilt/count_mat.H3K4me1.varfilt_1.5.rds"))
inf.k9me3 <- file.path(hubprefix, paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/varfilt/count_mat.H3K9me3.varfilt_1.5.rds"))

mat.k4me1 <- readRDS(inf.k4me1)
mat.k9me3 <- readRDS(inf.k9me3)

good.cells.k4me1 <- colnames(mat.k4me1)
good.cells.k9me3 <- colnames(mat.k9me3)


# Get good bins -----------------------------------------------------------

inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl.blfix/count_mat.H3K9me3-H3K4me1.filt_0.15_0.95_counts_and_l2r.blfix.rds"
mat.dbl <- readRDS(inf.dbl)

bins.keep <- 
bins.keep <- dat.bed$V4


# Write filterd mats ------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.single_match_dbl"
outf1 <- file.path(outdir, "count_mat.H3K4me1.match_dbl.rds")
outf2 <- file.path(outdir, "count_mat.H3K9me3.match_dbl.rds")

cols.keep.k4me1 <- colnames(mat.k4me1.raw) %in% good.cells.k4me1
cols.keep.k9me3 <- colnames(mat.k9me3.raw) %in% good.cells.k9me3

rows.keep.k4me1 <- rownames(mat.k4me1.raw) %in% bins.keep
rows.keep.k9me3 <- rownames(mat.k9me3.raw) %in% bins.keep

mat.k4me1.raw.filt <- mat.k4me1.raw[rows.keep.k4me1, cols.keep.k4me1]
mat.k9me3.raw.filt <- mat.k9me3.raw[rows.keep.k9me3, cols.keep.k9me3]

saveRDS(mat.k4me1.raw.filt, file = outf1)
saveRDS(mat.k9me3.raw.filt, file = outf2)




