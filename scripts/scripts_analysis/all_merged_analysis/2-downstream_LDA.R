# Jake Yeung
# Date of Creation: 2019-12-07
# File: ~/projects/scchic/scripts/scripts_analysis/all_merged_analysis/2-downstream_LDA.R
# Downstream LDA 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

# Load  -------------------------------------------------------------------

# all marks
jmark <- "H3K27me3"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.TRUE/ldaOut.B6BM_AllMerged_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = strsplit(cell, split = "-")[[1]][[5]])

# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

# get variance
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~experi)


# Find topics  ------------------------------------------------------------



