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

# 
# # Load singles from csv ---------------------------------------------------
# 
# indir.raw <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/countTablesAndRZr1only_ByChromo.NewFilters")
# assertthat::assert_that(dir.exists(indir.raw))
# infs.k4me1 <- list.files(indir.raw, pattern = ".*.K4me1.*binsize_50000.csv", full.names = TRUE)
# infs.k9me3 <- list.files(indir.raw, pattern = ".*.K9me3.*binsize_50000.csv", full.names = TRUE)
# 
# mat.k4me1.raw.lst <- ReadMatSlideWinFormat(infs.k4me1)
# mat.k9me3.raw.lst <- ReadMatSlideWinFormat(infs.k9me3)
# 
# # Load singles from filtered ------------------------------------------------------------
# 
# inf.k4me1 <- file.path(hubprefix, paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/varfilt/count_mat.H3K4me1.varfilt_1.5.rds"))
# inf.k9me3 <- file.path(hubprefix, paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/varfilt/count_mat.H3K9me3.varfilt_1.5.rds"))
# 
# mat.k4me1 <- readRDS(inf.k4me1)
# mat.k9me3 <- readRDS(inf.k9me3)
# 
# good.cells <- c(colnames(mat.k4me1), colnames(mat.k9me3))

# Load doubles ------------------------------------------------------------

inf.dbl <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl.blfix/count_mat.H3K9me3-H3K4me1.filt_0.15_0.95_counts_and_l2r.blfix.rds"))
mat.dbl <- readRDS(inf.dbl)

good.bins <- rownames(mat.dbl)

good.bins.dat <- data.frame(bin = good.bins, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(Chromo = GetChromo(bin),
         Start = GetStart(bin),
         End = GetEnd(bin))



good.bins.bed <- subset(good.bins.dat, select = c(Chromo, Start, End, bin))

# Write good bins to output -----------------------------------------------

outf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl.blfix/good_bins_dbl.blfix.bed")

fwrite(good.bins.bed, file = outf, col.names = FALSE, sep = "\t")




