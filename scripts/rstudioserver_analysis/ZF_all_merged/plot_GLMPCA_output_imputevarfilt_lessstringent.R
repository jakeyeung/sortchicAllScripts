# Jake Yeung
# Date of Creation: 2020-04-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/plot_GLMPCA_output_imputevarfilt_lessstringent.R
#


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

# Load DE genes -----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/from_data/zebrafish.poisson.2020-04-27/diff_exprs_Chloe_seurat.full.ctypefilt.rds"
dat.de <- readRDS(inf.de)


# Contants ----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks
jmark <- "H3K4me1"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/GLMPCA_UMAP"

dir.create(outdir)

for (jmark in jmarks){
  
  print(jmark)
  
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster != "Unknown")  # remove the small cluster Unknown
  
  annot.glmpca.filt <- left_join(annot.glmpca.filt, subset(annot.louv, select = c(cell, var.imputed)))
  
  # plot uMAP 
  
  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m <- ggplot(annot.glmpca.filt, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + ggtitle("UMAP of LDA ")
  print(m)
  
  
  
  # Load GLMPCA -------------------------------------------------------------
  
  inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged.imputevarfilt.lessstringent/ZF_", jmark, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.GLMPCA_var_correction.mergebinsize_1000.binskeep_500.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.5.winsorize_TRUE.2020-04-29.RData")
  load(inf.glmpca, v=T)
  # do it on the GLMPCA? 
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  dat.umap.glm <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings) %>%
    left_join(., subset(annot.glmpca.filt, select = c(cell, cluster)))
  
  m <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(paste(jmark, "GLMPCA UMAP"))
  print(m)
  
  
  outpdf <- file.path(outdir, paste0(jmark, ".GLMPCA_output_UMAP.pdf")) 
  pdf(outpdf, useDingbats = FALSE)
  print(m)
  dev.off()
}
