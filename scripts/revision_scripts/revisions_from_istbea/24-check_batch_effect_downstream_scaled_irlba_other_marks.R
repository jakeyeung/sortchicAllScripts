# Jake Yeung
# Date of Creation: 2022-05-03
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/24-check_batch_effect_downstream_scaled_irlba_other_marks.R
# 

rm(list=ls())

library(irlba)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 60
jsettings[["min_dist"]] <- 0.5
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 4


jmark <- "k27me3"

jmarks <- c("k4me1", "k4me3", "k9me3", "k27me3"); names(jmarks) <- jmarks

ntopics <- 30
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"
for (jmark in jmarks){
  
# out <- parallel::mclapply(jmarks, function(jmark){
  print(jmark)
  
  # Load mat ----------------------------------------------------------------
  
  # inf.mat <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_batch_corrected.2022-04-19.rds"
  inf.mat <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_", jmark, "_batch_corrected.2022-05-02.rds")
  mat <- readRDS(inf.mat)
  
  # SVD ---------------------------------------------------------------------
  
  
  logodds <- mat  # alreay in -Inf Inf
  # remove mean and SVD
  logodds.centered <- t(scale(t(logodds), center = TRUE, scale = TRUE))
  # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
  # logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
  logodds.pca <- irlba(A = t(logodds.centered), nv = ntopics, scale = FALSE, center = FALSE)
  rownames(logodds.pca$u) <- colnames(logodds)
  rownames(logodds.pca$v) <- rownames(logodds)
  
  outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"
  outpdf <- file.path(outdir, paste0("plots_pca_output_mat_adj_irlba_scaled2.", jmark, ".", Sys.Date(), ".pdf"))
  saveRDS(logodds.pca, file = file.path(outdir, paste0("pca_output_mat_adj_irlba_scaled2.", jmark, ".", Sys.Date(), ".rds")))
  
  pdf(outpdf, useDingbats = FALSE)
  # logodds.pca <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/pca_output_mat_adj.2022-04-20.rds")
  
  PoV <- signif(logodds.pca$d^2/sum(logodds.pca$d^2), digits = 2)
  # summary(logodds.pca)
  
  U.init <- logodds.pca$u  # cells by k
  V.init <- logodds.pca$v  # genes by k, no need to transpose
  
  rownames(U.init) <- colnames(logodds)
  rownames(V.init) <- rownames(logodds)
  
  # saveRDS(logodds.pca, file = paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/pca_output_mat_adj.", Sys.Date(), ".rds"))
  
  # UMAP --------------------------------------------------------------------
  
  
  inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt"))
  dat.meta <- fread(inf.meta) %>%
    dplyr::select(cell, ctype.from.LL, batch, colcode, ctype) %>%
    rowwise() %>%
    mutate(cluster = ctype.from.LL,
           jrep2 = ifelse(batch == "New", "anew", "old"),
           cluster = ifelse(cluster == "HSCs", "aHSCs", cluster))
  
  
  dat.meta.colors <- subset(dat.meta, select = c(colcode, ctype.from.LL))
  dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]
  
  dat.umap.annot <- scchicFuncs::DoUmapAndLouvain(U.init, jsettings) %>%
    left_join(., dat.meta)
  
  # dat.umap.annot <- left_join(dat.umap, dat.meta.colors)
  
  m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(~batch) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    facet_wrap(~batch) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  dat.pca <- data.frame(cell = rownames(U.init), pc1 = U.init[, 1], pc2 =U.init[, 2], stringsAsFactors = FALSE) %>%
    left_join(., dat.meta) %>%
    rowwise() %>%
    mutate(ctype.order = ifelse(ctype.from.LL == "HSCs", "zzzHSCs", ctype.from.LL)) %>%
    mutate(ctype.order = ifelse(ctype.order == "LT", "zzLT", ctype.order)) %>%
    mutate(ctype.order = ifelse(ctype.order == "ST", "zST", ctype.order)) %>%
    arrange(ctype.order)
  
  m <- ggplot(dat.pca, aes(x = pc1, y = pc2, color = colcode)) + 
    geom_point() + 
    xlab(paste0("PC1 (", PoV[[1]], ")")) + 
    ylab(paste0("PC2 (", PoV[[2]], ")")) + 
    theme_bw() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(dat.pca, aes(x = pc1, y = pc2, color = colcode)) + 
    geom_point() + 
    xlab(paste0("PC1 (", PoV[[1]], ")")) + 
    ylab(paste0("PC2 (", PoV[[2]], ")")) + 
    theme_bw() + 
    facet_wrap(~batch) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  # plot all PCs
  
  jdims <- seq(ntopics - 1); names(jdims) <- jdims
  
  for (jdim1 in jdims){
    jdim2 <- jdim1 + 1
    dat.pca.tmp <- data.frame(cell = rownames(U.init), pc1 = U.init[, jdim1], pc2 =U.init[, jdim2], stringsAsFactors = FALSE) %>%
      left_join(., dat.meta) %>%
      rowwise() %>%
      mutate(ctype.order = ifelse(ctype.from.LL == "HSCs", "zzzHSCs", ctype.from.LL)) %>%
      mutate(ctype.order = ifelse(ctype.order == "LT", "zzLT", ctype.order)) %>%
      mutate(ctype.order = ifelse(ctype.order == "ST", "zST", ctype.order)) %>%
      arrange(ctype.order)
    
    m.tmp <- ggplot(dat.pca.tmp, aes(x = pc1, y = pc2, color = colcode)) + 
      geom_point() + 
      xlab(paste0("PC", jdim1, " (", PoV[[jdim1]], ")")) + 
      ylab(paste0("PC", jdim2, " (", PoV[[jdim2]], ")")) + 
      theme_bw() + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m.tmp)
    
    
    
    
  }
  
  
  
  dev.off()
  saveRDS(dat.pca, file = file.path(outdir, paste0("dat_pca_adj.", jmark, ".", Sys.Date(), ".rds")))
  saveRDS(dat.umap.annot, file = file.path(outdir, paste0("dat_umap_annot_adj.", jmark, ".", Sys.Date(), ".rds")))
  
  # return(list(dat.umap.annot = dat.umap.annot, dat.pca = dat.pca))
  
# }, mc.cores = length(jmarks))
} 
  
