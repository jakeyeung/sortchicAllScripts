# Jake Yeung
# Date of Creation: 2019-06-19
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/prepare_LDAs_for_projection.R
# Get LDA objects from previous analyses and save them to do projection later

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(umap)

kchoose <- 50
outdir <- "/Users/yeung/data/scchic/LDA_for_training"

# Give up: use the bamlist for merging directory directly 

# 
# # H3K4me3 -----------------------------------------------------------------
# 
# inf.h3k4me3 <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.stringent_filter/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50.Robj"
# load(inf.h3k4me3, v=T)
# 
# outname <- paste0(strsplit(basename(inf.h3k4me3), "\\.")[[1]][[2]], ".RData")
# 
# kvec <- lapply(out.lda, function(x) x@k)
# 
# indx <- which(kvec == kchoose)
# 
# lda.orig <- out.lda[[indx]]
# 
# save(lda.orig, file = file.path(outdir, outname))
# 
# 
# # H3K27me3 ----------------------------------------------------------------
# 
# 
# 
# 
# 
# 
# # H3K9me3 -----------------------------------------------------------------
# 
# 
# 
# # H3K4me1 -----------------------------------------------------------------
# 
# 
# jmark <- "H3K4me1"
# jbin <- TRUE
# kstr <- "25_30_35_50"
# inf.h3k4me1 <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.B6_H3K4me1_pcutoff_0.CountThres0.K-25_30_40_50.Robj"
# assertthat::assert_that(file.exists(inf.h3k4me1))
# load(inf.h3k4me1, v=T)
# (outname <- paste0(strsplit(basename(inf.h3k4me1), "\\.")[[1]][[2]], ".RData"))
# 
# kvec <- lapply(out.lda, function(x) x@k)
# indx <- which(kvec == kchoose)
# lda.orig <- out.lda[[indx]]
# 
# save(lda.orig, file = file.path(outdir, outname))
# 
# 
# 
