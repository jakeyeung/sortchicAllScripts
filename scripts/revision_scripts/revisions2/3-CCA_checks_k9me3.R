# Jake Yeung
# Date of Creation: 2022-07-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/3-CCA_checks_k9me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k4me1", "k9me3", "k4me3", "k27me3"); names(jmarks) <- jmarks
# Metas -------------------------------------------------------------------


inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)



dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

dat.meta.merge <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  subset(jdat, select = c(cell, ctype.from.LL, colcode)) %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

# Check CCA outliers?  ----------------------------------------------------

# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins/UMAP_of_CCA_repressive_k9me3_and_k4me1.2022-07-22.RData"
jmarkref <- "k9me3"
jmarkothers <- jmarks[jmarks != jmarkref]
dat.umap.cca.lst <- lapply(jmarkothers, function(jmarkother){
  print(jmarkother)
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500/UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".2022-07-22.RData")
  load(inf.meta, v=T)
  dat.umap.cca <- dat.umap.cca %>%
    left_join(., dat.meta.merge)
})

m.lst <- lapply(jmarkothers, function(jmarkother){
  ggplot(dat.umap.cca.lst[[jmarkother]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    ggtitle(paste(jmarkref, "in", jmarkother)) + 
    scale_color_identity() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst)

print(m.lst$k4me1)
bad.cells.k4me1 <- subset(dat.umap.cca.lst$k4me1, umap1 < -5 | umap2 > 20)$cell

outbadcells <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/k4me1_bad_cells.txt"
readr::write_lines(bad.cells.k4me1, file = outbadcells)

ggplot(dat.umap.cca.lst$k4me1 %>% filter(!cell %in% bad.cells.k4me1), aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  scale_color_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  
# Can we get UMAP from k9me3 LDA ------------------------------------------

# inf.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins/lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins.H3K9me3.2022-07-21/ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins.H3K9me3.2022-07-21.Robj"
inf.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins_keeptop_500/lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_500.H3K9me3.2022-07-22/ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_500.H3K9me3.2022-07-22.Robj"
load(inf.k9, v=T)

tm <- out.lda %>%
  posterior() %>%
  scchicFuncs::AddTopicToTmResult()


library(hash)
library(igraph)
library(umap)
jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

dat.umap.pca <- DoUmapAndLouvain(tm$topics, jsettings)
dat.umap.pca.annot <- dat.umap.pca %>%
  left_join(., dat.meta.merge)
  
ggplot(dat.umap.pca.annot, aes(x = umap1, y = umap2, color = colcode))  + 
  geom_point() + 
  theme_bw() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

