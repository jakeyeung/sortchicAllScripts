# Jake Yeung
# Date of Creation: 2022-07-26
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/7-summarize_CCA_plots_all_pairs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/primetime_plots"

outpdf <- file.path(outdir, paste0("CCA_summary_all_pairs_k9me3_no_TSS.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks

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


# Plot the easy ones ------------------------------------------------------

# # actives
# inf.actives <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/CCA_k4me1_and_k4me3.repressive.2022-07-21.RData"
# load(inf.actives, v=T)
# 
# 

indir.cca <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/refmarks_merged"
jpairs <- list("k4me1_and_k4me3" = c("k4me1", "k4me3"), 
              "k4me1_and_k27me3" = c("k4me1", "k27me3"), 
              "k4me3_and_k27me3" = c("k4me3", "k27me3"))

dat.umap.cca.lst <- lapply(jpairs, function(jpair){
  jmarkref <- jpair[[1]]
  jmarkother <- jpair[[2]]
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500_bymark_factor_-1/UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".2022-07-22.RData")
  # inf.meta <- file.path(indir.cca, paste0("UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".2022-07-21.RData"))
  inf.meta <- file.path(indir.cca, paste0("UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".RData"))
  load(inf.meta, v=T)
  dat.umap.cca <- dat.umap.cca %>%
    left_join(., dat.meta.merge)
})


# Plot pairs --------------------------------------------------------------

jpairs.names <- names(jpairs); names(jpairs.names) <- jpairs.names

m.lst.easy <- lapply(jpairs.names, function(jpair.name){
  ggplot(dat.umap.cca.lst[[jpair.name]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jpair.name) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})
m.lst.easy


m.lst.easy.bymark <- lapply(jpairs.names, function(jpair.name){
  ggplot(dat.umap.cca.lst[[jpair.name]], aes(x = umap1, y = umap2, color = mark)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(jpair.name) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})
m.lst.easy.bymark



# Show k9me3  -------------------------------------------------------------

indir.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500_bymark_factor_-1"
# indir.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500_bymark_factor_1_with_TSS"
jmarkref <- "k9me3"
jmarkothers <- jmarks[jmarks != jmarkref]

dat.umap.cca.k9.lst <- lapply(jmarkothers, function(jmarkother){
  print(jmarkother)
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_500_bymark_factor_-1/UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkother, ".2022-07-22.RData")
  load(inf.meta, v=T)
  dat.umap.cca <- dat.umap.cca %>%
    left_join(., dat.meta.merge)
})

m.k9.lst <- lapply(jmarkothers, function(jmarkother){
  m <- ggplot(dat.umap.cca.k9.lst[[jmarkother]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jmarkother) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})


m.k9.bymark.lst <- lapply(jmarkothers, function(jmarkother){
  m <- ggplot(dat.umap.cca.k9.lst[[jmarkother]], aes(x = umap1, y = umap2, color = mark)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(jmarkother) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})


# Plot some loadings ?  ---------------------------------------------------


# k4me1 and k27me3 
# jpair.check <- "k4me1_and_k4me3"
jpair.check <- "k4me1_and_k27me3"
for (jpair.check in jpairs.names){
  jinf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/refmarks_merged/CCA_", jpair.check, ".repressive.RData")
  assertthat::assert_that(file.exists(jinf))
  load(jinf, v=T)
  mat.cca <- marks.combined@reductions$cca@cell.embeddings
  
  # plot first two components
  dat.cca <- data.frame(cell = rownames(mat.cca), mat.cca, stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.merge)
  
  ctypes.ordered <- c("Granulocytes", "Bcells", "Eryths", "NKs", "Monocytes", "DCs", "Basophils", "MEP", "CMP", "pDCs", "MPPs", "LT", "ST", "HSCs")
  dat.cca$ctype.from.LL <- factor(dat.cca$ctype.from.LL, levels = ctypes.ordered)
  
  dat.cca <- dat.cca %>%
    arrange(ctype.from.LL)
  
  m <- ggplot(dat.cca, aes(x = CC_1, y = -1 * CC_2, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jpair.check) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m.split <- ggplot(dat.cca, aes(x = CC_1, y = -1 * CC_2, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jpair.check) + 
    theme_bw() + 
    facet_wrap(~ctype.from.LL) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.split)
  
  
  m2 <- ggplot(dat.cca, aes(x = CC_2, y = CC_3, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jpair.check) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m2)
  
    
  m3 <- ggplot(dat.cca, aes(x = CC_3, y = CC_4, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jpair.check) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m3)
  
}

print(m.k9.lst)

print(m.k9.bymark.lst)
dev.off()