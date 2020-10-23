# Jake Yeung
# Date of Creation: 2020-08-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/K562_analysis/compare_K562_projections.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Plot projections --------------------------------------------------------

# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/projections_K562_spikein_integrate_with_old.top5000/projections.count_mats_old_binsize_50000_genomewide.top5000.H3K4me3.old.K-30/count_mats_old_binsize_50000_genomewide.top5000.projection.RData"
inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/projections_K562_spikein_integrate_with_old/projections.count_mats_old_binsize_50000_genomewide.H3K27me3.old.K-30/count_mats_old_binsize_50000_genomewide.projection.RData"
load(inf, v=T)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# dat.impute.log.old <- t(tm.r)

umap.out <- umap(posterior(out.objs$out.lda)$topics, jsettings)
dat.umap <- DoUmapAndLouvain(posterior(out.objs$out.lda)$topics, jsettings = jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# add new 

topics.new <- out.lda.predict$topics

predout <- predict(umap.out, data = topics.new)
dat.pred <- data.frame(umap1 = predout[, 1], umap2 = predout[, 2], cell = rownames(predout), stringsAsFactors = FALSE)
dat.pred$experi <- "VAN4969"

dat.orig <- data.frame(cell = dat.umap$cell, umap1 = dat.umap$umap1, umap2 = dat.umap$umap2, experi = "orig", stringsAsFactors = FALSE)

dat.merge <- bind_rows(dat.orig, dat.pred)

ggplot(dat.pred, aes(x = umap1, y = umap2, color = experi)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = umap1, y = umap2, color = experi)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add UMAP  ---------------------------------------------------------------


