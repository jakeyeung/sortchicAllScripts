# Jake Yeung
# Date of Creation: 2019-06-03
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/analysis_checks/downsampled_variance_analysis.R
# Analyze downsampled variance over pseudotime
# 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(topicmodels)
library(umap)
library(Matrix)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/PlotFunctions.R")

# inf.down <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.downsamp_2000/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50.Robj"
# inf.down <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.downsamp_700/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50.Robj"
inf.down <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.stringent_filter/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50.Robj"
load(inf.down, v=T)

out.lda <- out.lda[[4]]


# Plot UMAP  --------------------------------------------------------------

tm.result <- posterior(out.lda)

out.umap <- umap(tm.result$topics)

dat.umap.long <- data.frame(umap1 = out.umap$layout[, 1], umap2 = out.umap$layout[, 2], cell = rownames(out.umap$layout))

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_classic() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Calculate variance across cells -----------------------------------------

jfac <- 10^6
jpseudo <- 0
dat.mat <-  t(tm.result$terms) %*% t(tm.result$topics)
# log2 transform
dat.mat <- log2(dat.mat * jfac + jpseudo)

cells.sd <- GetCellSd(dat.mat, "", log2.scale = FALSE, fn = var) 

# merge with umap and plot

dat.merged <- left_join(dat.umap.long, cells.sd) %>%
  dplyr::rename(cell.var = cell.sd)

PlotXYWithColor(dat.merged, xvar = "umap1", yvar = "umap2", cname = "cell.var", jsize = 5)

# plot(density(Matrix::colSums(count.mat)))