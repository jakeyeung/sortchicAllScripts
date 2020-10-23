# Jake Yeung
# Date of Creation: 2020-09-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/project_LDA_downstream.R
# description


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

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_projection_onto_old.VAN5046/ldaOut.BM_H3K4me1_varfilt_countmat.2020-02-11.AllMerged.K-30.x.count_mat_H3K4me1_l2r_filt.2020-09-13.minl2r_-1.varfilt_1.RData"
load(inf, v=T)


# Plot original  ----------------------------------------------------------

DoUmapFromLDA

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

tm.result.orig <- posterior(out.objs$out.lda)
topics.mat <- tm.result.orig$topics
dat.umap <- umap(topics.mat, jsettings)

rownames(dat.umap$layout) <- rownames(topics.mat)
dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long$experi <- "old"

# Add projections ---------------------------------------------------------


dat.umap.pred <- predict(dat.umap, data = out.lda.predict$topics)

dat.umap.long.pred <- as.data.frame(dat.umap.pred)
colnames(dat.umap.long.pred) <- c("umap1", "umap2")
dat.umap.long.pred$cell <- rownames(dat.umap.long.pred)
dat.umap.long.pred$experi <- "new"

dat.merge <- rbind(dat.umap.long, dat.umap.long.pred)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() + 
  facet_wrap(~experi) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


