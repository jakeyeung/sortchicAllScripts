# Jake Yeung
# Date of Creation: 2019-11-28
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_mergetagbams.R
# Merge bams before tagging to see difference in LDA

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

library(irlba)

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged/lda_outputs.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-27.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-27.K-30.Robj"
load(inf, v=T)


tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


lsi.out <- RunLSI(as.matrix(BinarizeMatrix(count.mat)))
umap.out.lsi <- umap(lsi.out$u, config = jsettings)
dat.umap.long.lsi <- data.frame(cell = rownames(umap.out.lsi$layout), umap1 = umap.out.lsi$layout[, 1], umap2 = umap.out.lsi$layout[, 2], stringsAsFactors = FALSE)
ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


