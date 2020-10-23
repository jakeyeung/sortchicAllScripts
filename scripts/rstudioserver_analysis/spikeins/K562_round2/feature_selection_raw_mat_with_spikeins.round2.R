# Jake Yeung
# Date of Creation: 2020-09-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/feature_selection_raw_mat_with_spikeins.round2.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)

library(scchicFuncs)

topn <- 5000

# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/merged"
outdir <- file.path(indir, "genes_filt")
dir.create(outdir)

infs <- list.files(indir, pattern = "K562.*.rds", all.files = TRUE, full.names = TRUE)

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.RData"
load(inf.spikeins, v=T)

lapply(infs, function(inf){
  bname <- ClipLast(basename(inf), jsep = ".rds")
  outname <- paste0(bname, ".topn_", topn, ".rds")
  outnamepdf <- paste0(bname, ".topn_", topn, ".pdf")
  outf <- file.path(outdir, outname)
  outpdf <- file.path(outdir, outnamepdf)
  assertthat::assert_that(!file.exists(outf))
  assertthat::assert_that(!file.exists(outpdf))
  
  
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
})

