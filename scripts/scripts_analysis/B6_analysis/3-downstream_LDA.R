# Jake Yeung
# Date of Creation: 2019-05-08
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/3-downstream_LDA.R
# Check out LDA outputs

rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load files --------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
  
  jmark <- "H3K4me1"
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/lda_out_meanfilt.B6_", jmark, "_binfilt_cellfilt.CountThres0.K-25_30_35_50.Robj")
  assertthat::assert_that(file.exists(inf))
  
  
  # Process LDA -------------------------------------------------------------
  
  out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE)
  # load(inf, v=T)
  # out.lda <- 
  
  # Plot UMAP ---------------------------------------------------------------
  
  jmetric.louv='euclidean'
  jmindist.louv=0.3
  jseed.louv=123
  
  nn.louv.new <- c(28, 35, 33, 31)
  jmindist.new <- c(0.2, 0.2, 0.3, 0.3)
  nn.new <- c(40, 48, 45, 27)
  custom.settings.new.lst <- mapply(function(x, y) GetUmapSettings(x, jmetric.louv, y, jseed.louv), nn.new, jmindist.new, SIMPLIFY = FALSE)
  custom.settings.louv.new.lst <- lapply(nn.louv.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))
  
  names(custom.settings.new.lst) <- jmarks
  names(custom.settings.louv.new.lst) <- jmarks
  
  dat.umap <- umap(out.objs$tm.result$topics, config = custom.settings.new.lst[[jmark]])
  
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2])
  
  # PlotXYWithColor(dat.umap.long, xvar = "umap1", yvar = "umap2")
  
  ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Find interesting topics -------------------------------------------------


