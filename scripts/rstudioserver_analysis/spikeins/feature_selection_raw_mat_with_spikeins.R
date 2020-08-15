# Jake Yeung
# Date of Creation: 2020-08-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/feature_selection_raw_mat_with_spikeins.R
# Do a feature selection using spikeins as "total" 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)

library(scchicFuncs)

topn <- 5000
# jsuffix <- "G1_G2_S."
# jsuffix <- "G1only."
jsuffix <- ""

# Load data  --------------------------------------------------------------

# jmark <- "H3K4me1"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, ".", jsuffix, "rds")
  outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, ".", jsuffix, "topn_", topn, ".rds")
  outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, ".", jsuffix, "topn_", topn, ".pdf")
  assertthat::assert_that(!file.exists(outf))
  
  
  inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData"
  
  load(inf.spikeins, v=T)
  
  mat <- readRDS(inf)
  
  # calculate devaicne for ever yrow
  
  dat.spikeins.sub <- dat.spikeins.mat[colnames(mat), ]
  xvec <- mat[1, ]
  p <- sum(xvec) / sum(dat.spikeins.sub$spikeincounts)
  nvec <- dat.spikeins.sub$spikeincounts
  
  scchicFuncs::binomial_deviance(x = xvec, p = p, n = nvec)
  
  
  
  # Rank genes by deviance --------------------------------------------------
  
  gdevs <- apply(mat, 1, function(xvec){
    scchicFuncs::binomial_deviance(x = xvec, p = sum(xvec) / sum(nvec), n = nvec)
  })
  
  pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
  plot(density(gdevs), main = "All genes", xlab = "Binomial Deviance")
  
  # check hihg dev gene
  topg.i <- which.max(gdevs)
  topg <- rownames(mat)[topg.i]
  cdat <- data.frame(count = mat[topg, ], cell = colnames(mat), stringsAsFactors = FALSE)
  
  plot(density(cdat$count), main = topg, xlab = "counts")
  
  dev.off()
  
  # Take top 5000 genes? ----------------------------------------------------
  
  genes.keep.vec <- sort(gdevs, decreasing = TRUE)[1:topn]
  genes.keep <- names(genes.keep.vec)
  
  mat.filt <- mat[genes.keep, ]
  
  
  # Write to output  --------------------------------------------------------
  print("Final dimension")
  print(dim(mat.filt))
  
  saveRDS(mat.filt, file = outf)
}

