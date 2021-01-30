# Jake Yeung
# Date of Creation: 2021-01-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/10-check_LDA_k9dynamicbins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettingsspread <- 8

hubprefix <- "/home/jyeung/hub_oudenaarden"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load metas --------------------------------------------------------------

indir.metas <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf.metas <- file.path(indir.metas, fname)
  fread(inf.metas)
})


# Load LDA  ---------------------------------------------------------------


dat.umaps.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  
  inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_dynamic_bins/lda_outputs.count_mat_k9_dynamic_bins_50kb.", jmark, ".2021-01-28.K-30.binarize.FALSE/ldaOut.count_mat_k9_dynamic_bins_50kb.", jmark, ".2021-01-28.K-30.Robj"))
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings) %>%
    left_join(., subset(dat.metas[[jmark]], select = c(cell, cluster, clustercol)))
  return(dat.umap)
})

jmark.tmp <- "H3K4me1"
jmark.tmp <- "H3K4me3"
jmark.tmp <- "H3K27me3"

jmark.tmp <- "H3K9me3"
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umaps.lst[[jmark.tmp]], aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



