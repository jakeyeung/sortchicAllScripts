# Jake Yeung
# Date of Creation: 2020-11-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/10-check_glmpca_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(JFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load glmpca  -------------------------------------------------------------

# jmark <- "H3K27me3"
jmark <- "H3K9me3"
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# 
# indir.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs"
# # inf.glmpca <- file.path(indir.glmpca, paste0("glmpca_", jmark, ".RData"))
# inf.glmpca <- file.path(indir.glmpca, paste0("glmpca.", jmark, ".bincutoff_10000.binskeep_0.RData"))
# 
# load(inf.glmpca, v=T)
# 
# dat.glmpca <- DoUmapAndLouvain(glm.out$factors, jsettings)
# 
# # dat.glmpca <- data.frame(cell = rownames(glm.out$factors), )
# # dat.glmpca <- DoLouvain(topics.mat = glm.out$factors, custom.settings.louv = jsettings, dat.umap.long = )

# 
# ggplot(dat.glmpca, mapping = aes(x = umap1, y = umap2, color = louvain)) + 
#   geom_point() + 
#   scale_color_manual(values = cbPalette) + 
#   theme_bw() + 
#   ggtitle(jmark, "GLMPCA") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Compare with LDA ?  -----------------------------------------------------

indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins"
fname.lda <- paste0("lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj")
inf.lda <- file.path(indir.lda, fname.lda)

load(inf.lda, v=T)

tm.result <- posterior(out.lda)
dat.lda <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

ggplot(dat.lda, mapping = aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  ggtitle(jmark, "LDA") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add celltype info -------------------------------------------------------

indir.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04"
fname.annot <- paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt")

dat.annot <- fread(file.path(indir.annot, fname.annot))

# dat.glmpca.annot <- left_join(dat.glmpca, subset(dat.annot, select = c(cell, cluster)))
dat.lda.annot <- left_join(dat.lda, subset(dat.annot, select = c(cell, cluster)))

# m.glmpca <- ggplot(dat.glmpca.annot, aes(x = umap1, y = umap2, color = cluster)) + 
#   geom_point() +  
#   scale_color_manual(values = cbPalette) + 
#   ggtitle(jmark, "GLMPCA") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.lda <- ggplot(dat.lda.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() +  
  scale_color_manual(values = cbPalette) + 
  ggtitle(jmark, "LDA") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

print(m.lda)

# multiplot(m.glmpca, m.lda, cols = 2)


# Do stuff  ---------------------------------------------------------------


