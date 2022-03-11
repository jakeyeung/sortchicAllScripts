# Jake Yeung
# Date of Creation: 2022-01-29
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/4-glmpca_downstream_varfilt2.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)


library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

# jmark <- "k27me3"
jiter <- "500"
jiter <- "1000"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/glmpca_outputs/varfilt/iters"
dir.create(outdir)
# jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarks <- c("k9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  print(jmark)
  
  # Load glmpca  ------------------------------------------------------------
  
  # inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_500.RData")
  inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_", jiter, ".RData")
  load(inf.glmpca, v=T)
  
  bname <- basename(inf.glmpca)
  pdfname <- paste0(ClipLast(x = bname, jsep = "\\.", jsep.out = "."), ".pdf")
  outpdf <- file.path(outdir, pdfname)
  
  
  # Load meta  --------------------------------------------------------------
  
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/metadata_plate_experi_batch._dynamic_bins_merged_with_old.", jmark, ".2022-01-28.txt")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  dat.meta <- fread(inf.meta)
  
  
  # Show UMAP  --------------------------------------------------------------
  
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  
  dat.meta.filt <- subset(dat.meta, select = c(cell, ctype, batch, Batch, plate))
  dat.umap.annot <- left_join(dat.umap, dat.meta.filt)
  
  m.umap <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m.umap2 <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  # Show individual glmpcas -------------------------------------------------
  
  glm.factors <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.filt)
  
  dims <- seq(ncol(glm.out$factors) - 1)
  # jfac1 <- "dim3"
  # jfac2 <- "dim4"
  
  pdf(outpdf, useDingbats = FALSE)
  print(m.umap)
  print(m.umap2)
  for (d1 in dims){
    d2 <- d1 + 1
    jfac1 <- paste0("dim", d1)
    jfac2 <- paste0("dim", d2)
    
    m1 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "ctype")) + 
      geom_point() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
    m2 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "ctype")) + 
      geom_point() + 
      theme_bw() + 
      facet_wrap(~ctype) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    
    print(m1)
    print(m2)
  }
  dev.off()
  
  
  
}
