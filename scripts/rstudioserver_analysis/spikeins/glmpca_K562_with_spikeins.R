# Jake Yeung
# Date of Creation: 2020-08-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/glmpca_K562_with_spikeins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load spikeincounts ------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
infs.chromo <- list.files(indir, pattern = "K562-EtOH-.*.csv", full.names = TRUE)

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos <- lapply(infs.chromo, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = NA)
}) %>%
  bind_rows() %>%
  filter(chromo == "1")

dat.spikeins.mat <- as.data.frame(subset(dat.chromos, select = c(samp, spikeincounts)))
rownames(dat.spikeins.mat) <- dat.spikeins.mat$samp

# Load data ---------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562"
mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inname <- paste0("K562_count_tables_50000.", jmark, ".G1_G2_S.rds")
  inf <- file.path(indir, inname)
  readRDS(inf)
})


Y <- as.matrix(mats.lst$H3K4me3)
spikeincounts.sub <- dat.spikeins.mat[colnames(Y), ]

spikeincounts <- spikeincounts.sub$spikeincounts; names(spikeincounts) <- spikeincounts.sub$samp


assertthat::assert_that(ncol(Y) == length(spikeincounts))

system.time(
  glmpcaout <- glmpca::glmpca(Y = Y, L = 30, fam = "poi", penalty = 1, sz = spikeincounts)
)

save(glmpcaout, Y, spikeincounts, file = outf)


