# Jake Yeung
# Date of Creation: 2020-10-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/check_fit_inputs.R
# Check i removed properly the clusters

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarkdbl <- "H3K4me1xH3K9me3"
remove.na <- TRUE

# outsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters"
outsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters"

indir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_", outsuffix)
assertthat::assert_that(dir.exists(indir))

infobjs <- file.path(indir, paste0("SetupObjs_AllMerged_", outsuffix, ".clstr_by_louvain_", jmarkdbl, ".removeNA_", remove.na, ".RData"))
load(infobjs, v=T)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")

ggplot(dat.louv$H3K4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.louv$H3K9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

