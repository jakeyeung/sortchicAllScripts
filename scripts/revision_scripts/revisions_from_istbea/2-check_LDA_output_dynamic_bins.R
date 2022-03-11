# Jake Yeung
# Date of Creation: 2022-01-27
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/2-check_LDA_output_dynamic_bins.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(hash)
library(igraph)
library(umap)
library(scchicFuncs)


jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123


jmark <- "k4me1"

jsuffix1 <- "_dynamic_bins"
jsuffix2 <- "_merged_with_old"

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata"

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3")

for (jmark in jmarks){
  print(jmark)
  outf.meta <- file.path(outdir, paste0("metadata_plate_experi_batch.", jsuffix1, jsuffix2, ".", jmark, ".", Sys.Date(), ".txt"))
  outf.pdf <- file.path(outdir, paste0("LDA_outputs", jsuffix1, jsuffix2, ".", jmark, ".", Sys.Date(), ".pdf"))
  indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed", jsuffix1)
  inf <- file.path(indir, paste0("ldaAnalysis_fripfilt_BM_", jmark, "/lda_outputs.count_mat", jsuffix2, ".", jmark, ".2022-01-26/ldaOut.count_mat", jsuffix2, ".", jmark, ".2022-01-26.Robj"))
  assertthat::assert_that(file.exists(inf))
  
  pdf(outf.pdf, useDingbats = FALSE)
  
  load(inf, v=T)
  
  tm <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm$topics, jsettings = jsettings)
  
  dat.umap <- dat.umap %>%
    rowwise() %>%
    mutate(plate = ClipLast(x = cell, jsep = "_"), 
           experi = ClipLast(x = plate, jsep = "-"), 
           batch = ifelse(grepl("PZ-sortChIC-BM-SL", experi), "New", "Old"), 
           mark = jmark)
  
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = plate)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~experi) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  # save output metafiles
  fwrite(dat.umap, file = outf.meta)
  
  dev.off()
  
  
  
  
}
