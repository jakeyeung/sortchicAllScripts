# Jake Yeung
# Date of Creation: 2020-10-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.celltyping.R
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
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}



# Functions ---------------------------------------------------------------




# Constnats ---------------------------------------------------------------


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")

# jmarks <- c("H3K27me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load presaved obj -------------------------------------------------------



inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all.merge_with_old/LDA_downstream_objects.2020-10-18.again.RData"

load(inf, v=T)



# Check plates ------------------------------------------------------------


jmark <- "H3K27me3"

jsub <- subset(dat.merge, mark == jmark)


ggplot(jsub, aes(x = umap1, y = umap2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 

ggplot(jsub %>% filter(grepl("rep3-H3K27me3-1$", plate)), aes(x = umap1, y = umap2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 

