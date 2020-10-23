# Jake Yeung
# Date of Creation: 2020-09-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/2-check_mouse_BM_glmpca_H3K4me1_H3K4me3_H3K27me3_together.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(glmpca)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette2 <- cbPalette
cbPalette2[[8]] <- cbPalette[[5]]
cbPalette2[[5]] <- cbPalette[[8]]

# Load annots BM ----------------------------------------------------------

# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
jmark <- "H3K4me1"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# jmarks

dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  indir.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/GLMPCA_outputs.Faster.pois.RemoveSmallPeaksKeepBestPlates2.LDAfromPeaks")
  dir.exists(indir.glmpca)
  inf.glmpca <- file.path(indir.glmpca, paste0("PZ_", jmark, ".AllMerged.KeepBestPlates2.LDAfromPeaks.GLMPCA_var_correction.mergebinsize_100.binskeep_250.covar_cell.var.within.sum.norm.log2.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-08-26.RData"))
  assertthat::assert_that(file.exists(inf.glmpca))
  
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables/BM_AllMerged.", jmark, ".cell_cluster_table.txt")
  
  load(inf.glmpca, v=T)
  
  dat.var <- subset(dat.merge2, select = c(cell, cell.var.within.sum.norm))
  
  # Load GLMPCA output ------------------------------------------------------
  
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings) %>%
    left_join(., dat.merge2, by = "cell")
  
  
  dat.annot <- fread(inf.annot)
  
  dat.umap.annot <- left_join(dat.umap, dat.annot, by = "cell")
  return(dat.umap.annot)
})

head(dat.umap.annot.lst$H3K4me1)
head(dat.umap.annot.lst$H3K4me3)
head(dat.umap.annot.lst$H3K27me3)

ggplot(dat.umap.annot.lst$H3K4me3, aes(x = umap1.x, y = umap2.x, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette2, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Clean up K4me3 clusters --------------------------------------------------

dat.umap.annot.lst.filt <- lapply(dat.umap.annot.lst, function(jdat){
  subset(jdat, cluster != "")
})

dat.umap.annot.lst.filt$H3K4me3$cluster.renamed <- sapply(dat.umap.annot.lst.filt$H3K4me3$cluster, function(x){
  xsplit1 <- strsplit(x, split = "_")[[1]][[1]]
  xsplit2 <- strsplit(xsplit1, split = "-")[[1]][[1]]
  xsplit2 <- gsub("ILC", "InnateLymph", xsplit2)
  xsplit2 <- gsub("Dendritic", "DC", xsplit2)
  xsplit2 <- gsub("pDC", "DC", xsplit2)
  xsplit2 <- gsub("Eryth", "Erythroblasts", xsplit2)
  return(xsplit2)
})

m0 <- ggplot(dat.umap.annot.lst.filt$H3K4me3, aes(x = umap1.x, y = umap2.x, color = cluster.renamed)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values = cbPalette2, na.value = "grey85") +
  xlab("") + ylab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", axis.ticks = element_blank(), axis.text = element_blank()) + 
  ggtitle("H3K4me3")


print(m0)

# Reannotate K27me3 clusters ----------------------------------------------

inf.reannot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K27me3.fromChIX/H3K27me3_cluster_renamed_from_chix_projections.txt"
dat.reannot <- fread(inf.reannot) %>%
  filter(orig != "")  %>%
  filter(!grepl("TooFewCells$", proj)) %>%
  filter(orig != "Eryth-Gfi1-_topic17") %>%
  rowwise() %>%
  mutate(proj = ifelse(orig %in% c("HSCs-Tead1-_topic9", "LinnegCore2_topic26"), "HSCs", proj)) %>%
  mutate(proj = ifelse(proj == "Eosinophils", "Basophils", proj))


clstr.rename <- hash::hash(dat.reannot$orig, dat.reannot$proj)
clstr.rename[["Eryth-Slc7a6-_topic1"]] <- "Erythroblasts"
clstr.rename[["Eryth-Sox6-_topic6"]] <- "Erythroblasts"
clstr.rename[["Eryth-Gfi1-_topic17"]] <- "Erythroblasts"
clstr.rename[["LinnegIsland2_topic29"]] <- "Basophils"
clstr.rename[["Bcells_topic16"]] <- "Bcells"

  
dat.umap.annot.lst.filt$H3K27me3$cluster.renamed <- sapply(dat.umap.annot.lst.filt$H3K27me3$cluster, function(x) AssignHash(x, clstr.rename))

m1 <- ggplot(dat.umap.annot.lst.filt$H3K27me3, aes(x = umap1.x, y = umap2.x, color = cluster.renamed)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values = cbPalette2, na.value = "grey85") +
  xlab("") + ylab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", axis.ticks = element_blank(), axis.text = element_blank()) + 
  ggtitle("H3K27me3")

ggplot(dat.umap.annot.lst.filt$H3K27me3, aes(x = umap1.x, y = umap2.x, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette2) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Clean up K4me1 ----------------------------------------------------------

dat.umap.annot.lst.filt$H3K4me1$cluster.renamed <- sapply(dat.umap.annot.lst.filt$H3K4me1$cluster, function(x){
  xsplit1 <- strsplit(x, split = "_")[[1]][[1]]
  xsplit2 <- strsplit(xsplit1, split = "-")[[1]][[1]]
  xsplit2 <- gsub("ILC", "InnateLymph", xsplit2)
  xsplit2 <- gsub("Dendritic", "DC", xsplit2)
  xsplit2 <- gsub("pDC", "DC", xsplit2)
  xsplit2 <- gsub("Eryth", "Erythroblasts", xsplit2)
  return(xsplit2)
})

m2 <- ggplot(dat.umap.annot.lst.filt$H3K4me1, aes(x = umap1.x, y = umap2.x, color = cluster.renamed)) + 
  geom_point() + 
  theme_minimal() + 
  scale_color_manual(values = cbPalette2, na.value = "grey85") +
  xlab("") + ylab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", axis.ticks = element_blank(), axis.text = element_blank()) + 
  ggtitle("H3K4me1")



outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/GLMPCA_peaks_primetime/H3K4me1_H3K4me3_H3K27me3_glmpca_peaks_primetime.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)
  print(m2)
  print(m1)
  print(m0)
  JFuncs::multiplot(m2, m1, cols = 2)
  JFuncs::multiplot(m2, m0, m1, cols = 3)
dev.off()


# Write output  -----------------------------------------------------------

for (jmark in jmarks){
  outtxt <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/GLMPCA_peaks_primetime/H3K4me1_H3K4me3_H3K27me3_glmpca_peaks_primetime.", Sys.Date(), ".", jmark, ".txt")
  print(jmark)
  outdat <- dat.umap.annot.lst.filt[[jmark]]
  fwrite(outdat, file = outtxt)
}
