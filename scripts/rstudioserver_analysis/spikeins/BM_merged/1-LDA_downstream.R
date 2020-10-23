# Jake Yeung
# Date of Creation: 2020-10-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.R
# 

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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load LDA  ---------------------------------------------------------------

# jmark <- "H3K4me3"
jmark <- "H3K9me3"
inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all/lda_outputs.count_mat.", jmark, ".filt_0.15_0.95_counts_and_l2r.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".filt_0.15_0.95_counts_and_l2r.K-30.Robj")
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  # scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#  Load GLMPCA  -----------------------------------------------------------


