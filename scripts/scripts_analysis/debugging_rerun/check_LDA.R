# Jake Yeung
# Date of Creation: 2019-12-20
# File: ~/projects/scchic/scripts/scripts_analysis/debugging_rerun/check_LDA.R
# Check LDA output

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

library(Rtsne)

# Load LDA  ---------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"

load(inf, v=T)


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# Do umap  ----------------------------------------------------------------

tm.result <- posterior(out.lda)

tsne.out <- Rtsne(tm.result$topics)

dat.tsne.long <- data.frame(cell = rownames(tm.result$topics), umap1 = tsne.out$Y[, 1], umap2 = tsne.out$Y[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell),
         experi = gsub("Bl6BMSC", "B6BMSC", experi), 
         plate = ClipLast(cell, jsep = "_"))

ggplot(dat.tsne.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_manual(values = cbPalette)

ggplot(dat.tsne.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_manual(values = cbPalette)  + 
  facet_wrap(~plate)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell),
         experi = gsub("Bl6BMSC", "B6BMSC", experi), 
         plate = ClipLast(cell, jsep = "_"))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_manual(values = cbPalette)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~plate)

# get var
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.tsne.long <- left_join(dat.tsne.long, dat.var)
dat.umap.long <- left_join(dat.umap.long, dat.var)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~plate)

ggplot(dat.tsne.long, aes(x = umap1, y = umap2, color = cell.var.within.sum)) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)

dat.cellsizes <- data.frame(cell = colnames(count.mat), cellsize = colSums(count.mat) / 5)

dat.tsne.long <- left_join(dat.tsne.long, dat.cellsizes)
dat.umap.long <- left_join(dat.umap.long, dat.cellsizes)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = log10(cellsize))) + geom_point(alpha = 0.5)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_color_viridis_c(direction = 1) + 
  facet_wrap(~plate)


