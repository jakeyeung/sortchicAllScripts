# Jake Yeung
# Date of Creation: 2022-01-28
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/4-glmpca_downstream.R
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
jsettings[["n_neighbors"]] <- 50
jsettings[["min_dist"]] <- 0.1
jsettings[["spread"]] <- 8 
jsettings[["random_state"]] <- 123

# jmark <- "k27me3"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/glmpca_outputs/k4me1_new_ctype_from_LL"
dir.create(outdir)
jmarks <- c("k4me1"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  
  
  # Load glmpca  ------------------------------------------------------------
  
  inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
  load(inf.glmpca, v=T)
  
  jfactors <- glm.out$factors
  l2norms <- ColNorms(jfactors)
  l2norms.frac <- l2norms / sum(l2norms)
  l2norms.frac.pretty <- signif(l2norms.frac, digits = 2) * 100
  
  bname <- basename(inf.glmpca)
  pdfname <- paste0(ClipLast(x = bname, jsep = "\\.", jsep.out = "."), ".", Sys.Date(), ".pdf")
  outpdf <- file.path(outdir, pdfname)
  
  
  # Load meta  --------------------------------------------------------------
  
  # # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/metadata_plate_experi_batch._dynamic_bins_merged_with_old.", jmark, ".2022-01-28.txt")
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  # dat.meta <- fread(inf.meta)
  
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".2022-04-15.txt")
  dat.meta <- fread(inf.meta)
  
  dat.meta.colors <- subset(dat.meta, select = c(ctype.from.LL, colcode))
  dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]
  
  # Show UMAP  --------------------------------------------------------------
  
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  
  dat.meta.filt <- subset(dat.meta, select = c(cell, ctype, batch, colcode, ctype.from.LL))
  dat.umap.annot <- left_join(dat.umap, dat.meta.filt)
  
  m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m.colors <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  # Show individual glmpcas -------------------------------------------------
  
  glm.factors <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.filt)
  
  dims <- seq(ncol(glm.out$factors) - 1)
  # jfac1 <- "dim3"
  # jfac2 <- "dim4"
  
  pdf(outpdf, useDingbats = FALSE)
  print(m)
  print(m.colors)
  for (d1 in dims){
    d2 <- d1 + 1
    jfac1 <- paste0("dim", d1)
    jfac2 <- paste0("dim", d2)
    
    m1 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "colcode")) + 
      geom_point() + 
      xlab(paste(jfac1, "(", l2norms.frac.pretty[[jfac1]], ")")) + 
      ylab(paste(jfac2, "(", l2norms.frac.pretty[[jfac2]], ")")) + 
      theme_bw() + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
    m2 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "colcode")) + 
      geom_point() + 
      theme_bw() + 
      xlab(paste(jfac1, "(", l2norms.frac.pretty[[jfac1]], ")")) + 
      ylab(paste(jfac2, "(", l2norms.frac.pretty[[jfac2]], ")")) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      facet_wrap(~ctype) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    
    m3 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "colcode")) + 
      geom_point() + 
      theme_bw() + 
      xlab(paste(jfac1, "(", l2norms.frac.pretty[[jfac1]], ")")) + 
      ylab(paste(jfac2, "(", l2norms.frac.pretty[[jfac2]], ")")) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      facet_wrap(~batch) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    
    print(m1)
    print(m2)
    print(m3)
    
  }
  dev.off()
  
  
  
}
