# Jake Yeung
# Date of Creation: 2020-03-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/Blood_with_intestines/2-LDA_downstream.R
# Blood only no intestines

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)


# Functions ---------------------------------------------------------------



# Load data ---------------------------------------------------------------

jmarklong <- "H3K4me1"
jmark <- "k4me1"

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3")
jmarkslong <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarkslong) <- jmarks

for (jmark in jmarks){
  print(jmark)
  jmarklong <- jmarkslong[[jmark]]
  
  inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-03-06.var_filt.UnenrichedAndAllMerged.KeepBestPlates2AddBlood/lda_outputs.countmat.Blood.", jmark, ".varfilt_0.3.2020-03-06.MergeWithBM.K-30.binarize.FALSE/ldaOut.countmat.Blood.", jmark, ".varfilt_0.3.2020-03-06.MergeWithBM.K-30.Robj")
  load(inf, v=T)
  
  tm.result <- posterior(out.lda)
  
  topics.mat <- tm.result$topics
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  dat.umap <- DoUmapAndLouvain(topics.mat, jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_")) %>%
    mutate(cond = GetCondFromSamp.blood(cell, mark = jmarklong))
  dat.umap$cond <- factor(dat.umap$cond, levels = c("Blood", "Unenriched", "Linneg", "StemCell"))
  
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  
  outpdf <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/AllMerged_withBlood/LDA_output_AllMerged_with_blood.", jmark, ".pdf")
  pdf(outpdf, useDingbats = FALSE)
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = plate)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +  
    scale_color_manual(values = cbPalette) + 
    facet_wrap(~cond) + ggtitle(paste("LDA output for", jmark))
  print(m)
  dev.off()
  
}






