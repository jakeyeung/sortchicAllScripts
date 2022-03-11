# Jake Yeung
# Date of Creation: 2022-01-27
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/2-check_LDA_output.R
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
library(scchicFuncs)


jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123


jmark <- "k9me3"
# inf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed/ldaAnalysis_fripfilt_BM_k4me1/lda_outputs.count_mat_new_only.k4me1.2022-01-26/ldaOut.count_mat_new_only.k4me1.2022-01-26.Robj"
inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed/ldaAnalysis_fripfilt_BM_", jmark, "/lda_outputs.count_mat_merged_with_old_fewer2.", jmark, ".2022-01-26/ldaOut.count_mat_merged_with_old_fewer2.", jmark, ".2022-01-26.Robj")
load(inf, v=T)

tm <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm$topics, jsettings = jsettings)

dat.umap <- dat.umap %>%
  rowwise() %>%
  mutate(plate = ClipLast(x = cell, jsep = "_"), 
         experi = ClipLast(x = plate, jsep = "-"))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.umap, aes(x = umap1, y = umap2, color = plate)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark) + 
  facet_wrap(~experi) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

