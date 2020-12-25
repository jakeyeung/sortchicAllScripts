# Jake Yeung
# Date of Creation: 2020-12-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_techrepmerged/6-correct_batch_H3K27me3_reseq.check_umap.R
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

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load umaps  -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep1rep2rep3reseq.peaks.varfilt/glmpca.H3K27me3.bincutoff_0.binskeep_0.platename_jrep.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData")
load(inf, v=T)

tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2.k4me1_and_k4me3_genes.K27me3only/matrix_imputed_mark_H3K27me3.newmetadata.pDCstoBasophils.txt")
dat.meta <- fread(inf.meta)

dat.umap.annot <- left_join(dat.umap, subset(dat.meta, select = -c(louvain, umap1, umap2)))

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~jrep) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
