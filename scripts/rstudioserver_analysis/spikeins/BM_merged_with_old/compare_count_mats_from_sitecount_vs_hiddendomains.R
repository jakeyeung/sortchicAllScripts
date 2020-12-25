# Jake Yeung
# Date of Creation: 2020-11-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/compare_count_mats_from_sitecount_vs_hiddendomains.R
# Chec rownames



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "H3K4me1"
inf1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMAllMerged2.from_peaks/filtNAcells_allbins/count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.rds"))
inf2 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat/filtNAcells_allbins/count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.rds"))

assertthat::assert_that(file.exists(inf1))
assertthat::assert_that(file.exists(inf2))


dat1 <- readRDS(inf1)
dat2 <- readRDS(inf2)

rownames(dat1) <- sapply(rownames(dat1), function(x) strsplit(x, ";")[[1]][[1]])
rownames(dat2) <- sapply(rownames(dat2), function(x) strsplit(x, ";")[[1]][[1]])

common.rows <- intersect(rownames(dat1), rownames(dat2))



# Check if it matches the mat files?  -------------------------------------





