# Jake Yeung
# Date of Creation: 2021-02-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K9me3_deeper/1-check_clustering_H3K9me3_no_hspcs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.H3K9me3.2021-01-02.txt"
dat.meta <- fread(inf.meta)

# Load data  --------------------------------------------------------------

inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins.RemoveHSPCs/lda_outputs.count_name.H3K9me3.k4_k9_dynamic_bins.RemoveHSPCs.2021-02-01.K-30.binarize.FALSE/ldaOut.count_name.H3K9me3.k4_k9_dynamic_bins.RemoveHSPCs.2021-02-01.K-30.Robj"
load(inf.lda, v=T)

tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics

dat.umap <- DoUmapAndLouvain(topics.mat, jsettings)

dat.umap.annot <- left_join(dat.umap, subset(dat.meta, select = c(cell, cluster, jrep, stype)))

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jsub <- subset(dat.umap.annot, umap2 > 5)

unique(sapply(jsub$cell, function(x) ClipLast(x = x, jsep = "_")))

