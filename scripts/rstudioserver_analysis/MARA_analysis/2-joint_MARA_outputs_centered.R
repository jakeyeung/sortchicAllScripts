# Jake Yeung
# Date of Creation: 2020-03-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/1-explore_MARA_outputs_centered_vs_uncentered.R
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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# Load GLMPCA  ------------------------------------------------------------


# load GLMPCA from bins 
jexperi <- "AllMerged"
mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

dat.merged.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  
  inf.glm <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
  assertthat::assert_that(file.exists(inf.glm))
  load(inf.glm, v=T)
  
  inf.annots <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annots))
  load(inf.annots, v=T)
  
  mdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
  assertthat::assert_that(dir.exists(mdir))
  
  
  dat.glmpca.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)
  
  dat.glmpca.umap <- dat.glmpca.umap %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), "_"))
  
  dat.glmpca.umap$cond <- sapply(dat.glmpca.umap$cell, GetCondFromSamp, mark = jmark)
  dat.glmpca.umap$cond <- factor(dat.glmpca.umap$cond, levels = c("Unenriched", "Linneg", "StemCell"))
  dat.glmpca.umap.annot <- left_join(dat.glmpca.umap, subset(dat.umap.glm.fillNAs, select = c(cell, cluster, topic.weight)))
  
  # Load MARA output  ------------------------------------------------------
  
  mara.out <- LoadMARA(mdir, make.cnames = FALSE)
  
  act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
  colnames(act.mat.clean) <- mara.out$act.mat$motif
  act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
    ungroup() %>%
    mutate(cell = gsub("\\.", "-", cell))
  
  dat.merge <- left_join(dat.glmpca.umap.annot, act.mat.clean.dat, by = "cell")
  
  return(dat.merge)

})