# Jake Yeung
# Date of Creation: 2021-01-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/17-LADA_checks_H3K4me1_H3K9me3_dbl.R
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

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks.withdbl <- c("H3K4me1", "H3K9me3", "H3K4me1-H3K9me3"); names(jmarks.withdbl) <- jmarks.withdbl


# Metas  ------------------------------------------------------------------

indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})

# LDA  --------------------------------------------------------------------

# jmark <- "H3K4me1"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins"
dat.umap.lst <- lapply(jmarks.withdbl, function(jmark){
  print(jmark)
  dname <- paste0("lda_outputs.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-29.K-30.binarize.FALSE/ldaOut.count_name.", jmark, ".k4_k9_dynamic_bins.2021-01-29.K-30.Robj")
  inf <- file.path(indir, dname)
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  topics.mat <- tm.result$topics
  dat.umap <- DoUmapAndLouvain(topics.mat, jsettings)
  if (jmark != "H3K4me1-H3K9me3"){
    dat.metas.tmp <- subset(dat.metas[[jmark]], select = c(-umap1, -umap2, -louvain))
    dat.umap.annot <- dat.umap %>%
      left_join(., dat.metas.tmp) %>%
      mutate(louvain = cluster)
  } else {
    dat.umap.annot <- dat.umap 
  }
  return(dat.umap.annot)
})



# Plot UMAP  --------------------------------------------------------------


m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[3[]])


