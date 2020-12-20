# Jake Yeung
# Date of Creation: 2020-10-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/calculate_foldchange_from_peaks.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(irlba)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/glmpca_spikeins_hoescht.K562_G1_G2/glmpca_spikeins_hoescht_summaries.K562_G1_G2.2020-10-27.PeakCounts.RData"
load(inf, v=T)


# Fit each peak -----------------------------------------------------------

mats.counts.filt$H3K4me1[1:5, 1:5]

jmark <- "H3K4me1"
jpeak <- "Peak_1"
xdat <- data.frame(mats.count.filt[[jmark]], ,stringsAsFactors = FALSE)
dat.meta <- subset(dat.glmpca.peak.lst[[jmark]])




