# Jake Yeung
# Date of Creation: 2020-08-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/check_LDA_from_peaks.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)


# Load LDA  ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap.lst <- lapply(jmarks, function(jmark){
  inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks/lda_outputs.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.binarize.FALSE/ldaOut.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.Robj")
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  
  m1 <- ggplot(dat.umap, aes(x = umap1, y = umap2)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m1)
  dat.umap$mark <- jmark
  return(dat.umap)
})


# Load annots -------------------------------------------------------------

dat.annot.lst <- lapply(jmarks, function(jmark){
  inf.annot <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt"))
  assertthat::assert_that(file.exists(inf.annot))
  dat.annot <- fread(inf.annot)
  dat.annot$mark <- jmark
  return(dat.annot)
})


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m.lst <- lapply(jmarks, function(jmark){
  dat.merge <- left_join(dat.umap.lst[[jmark]], subset(dat.annot.lst[[jmark]], select = c(cell, cluster)))
  m1 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~mark)
  print(m1)
  return(m1)
})


