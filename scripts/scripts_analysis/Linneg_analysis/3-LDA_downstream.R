# Jake Yeung
# Date of Creation: 2019-06-19
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/3-LDA_downstream.R
# Analyze LDA downstream

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(umap)

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.CountThres0.K-30_35_50_mindist_0.4_mindist_processed_lda.Rdata"
assertthat::assert_that(file.exists(inf))

load(inf, v=T)


# Plot output  ------------------------------------------------------------

dim(out.objs$tm.result$topics)

settings <- umap.defaults
settings$n_neighbors <- 58
settings$min_dist <- 0.4
settings$random_state <- 123
umap.out <- umap(out.objs$tm.result$topics, config = settings)

umap.long <- data.frame(umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], batch = sapply(rownames(umap.out$layout), function(x) strsplit(x, split = "-")[[1]][[1]]))

ggplot(umap.long, aes(x = umap1, y = umap2, color = batch)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

