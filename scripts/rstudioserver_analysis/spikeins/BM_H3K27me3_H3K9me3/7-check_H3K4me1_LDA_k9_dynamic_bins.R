# Jake Yeung
# Date of Creation: 2021-01-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/7-check_H3K4me1_LDA_k9_dynamic_bins.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8


hubprefix <- "/home/jyeung/hub_oudenaarden"

# Get emtas ---------------------------------------------------------------


ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned")
assertthat::assert_that(dir.exists(indir.meta))
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  print(inf)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})

dat.metas.reordered <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    jdat.reordered <- dat.metas[[jmark]] 
    jdat.reordered$cluster <- factor(jdat.reordered$cluster, levels = ctypes)
    jdat.reordered <- jdat.reordered %>%
      arrange(cluster)
  } else {
    jdat.reordered <- dat.metas[[jmark]]
  }
  return(jdat.reordered)
})



# Load LDA  ---------------------------------------------------------------


jmark <- "H3K4me1"
inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_dynamic_bins/lda_outputs.count_mat_by_peaks.H3K4me1.overlap_dynamic_k9_bins.K-30.binarize.FALSE/ldaOut.count_mat_by_peaks.", jmark, ".overlap_dynamic_k9_bins.K-30.Robj"))
load(inf.lda, v=T)


# Get outputs -------------------------------------------------------------

tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

dat.umap.annot <- left_join(dat.umap, subset(dat.metas.reordered[[jmark]], select = c(-umap1, -umap2, -louvain)))

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




