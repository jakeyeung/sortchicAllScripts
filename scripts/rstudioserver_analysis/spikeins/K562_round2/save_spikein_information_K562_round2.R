# Jake Yeung
# Date of Creation: 2020-09-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/save_spikein_information_K562_round2.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# load data  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir.chromo.g1 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
indir.chromo.cc <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/cellcycle_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")

assertthat::assert_that(dir.exists(indir.chromo.g1))
assertthat::assert_that(dir.exists(indir.chromo.cc))

infs.chromo.g1 <- list.files(indir.chromo.g1, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
infs.chromo.cc <- list.files(indir.chromo.cc, pattern = "K562-EtOH-.*.csv", full.names = TRUE)


jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos.g1 <- lapply(infs.chromo.g1, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = jchromos)
}) %>%
  bind_rows()

dat.chromos.cc <- lapply(infs.chromo.cc, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = jchromos)
}) %>%
  bind_rows()

dat.chromos <- rbind(dat.chromos.g1, dat.chromos.cc)

chromocounts <- subset(dat.chromos, chromo == "1", select = c(samp, chromocounts, spikeincounts))

dat.spikeins.mat <- as.data.frame(chromocounts)

# add spikeins
rownames(dat.spikeins.mat) <- dat.spikeins.mat$samp




# save to output

outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.RData"

# saveRDS(object = dat.spikeins.mat, file = outf)
save(dat.spikeins.mat, file = outf)

