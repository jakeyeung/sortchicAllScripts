# Jake Yeung
# Date of Creation: 2022-02-18
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/8-downstream_LL_celltyping_repressive_cleaned.R
# description
rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k9me3", "k27me3"); names(jmarks) <- jmarks
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_repressive_cleaned_Tcells"
jsuffix <- "dynamicbins"
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_repressive_cleaned_Tcells"

# jmark <- jmarks[[1]]

for (jmark in jmarks){
  
  outmeta <- file.path(outdir, paste0("metadata_reannotate_from_LLmat_", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  outpdf <- file.path(outdir, paste0("plots_reannotate_from_LLmat_", jsuffix, ".", jmark, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  
  inf <- file.path(indir, paste0("LLmat_dynamicbins_", jmark, ".2022-02-23.rds"))
  
  LL.all <- readRDS(inf)
  
  LL.mat <- t(LL.all)
  
  indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned"
  inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".dynamicbins.2022-02-18.txt"))
  dat.meta <- fread(inf.meta)
  
  ctypes.indx <- colnames(LL.mat)
  
  LL.best <- apply(LL.mat, 1, function(jrow){
    ctypes.indx[which.max(jrow)]
  })
  
  dat.LL.best <- data.frame(cell = names(LL.best), ctype.from.LL = LL.best, stringsAsFactors = FALSE)
  
  dat.meta.reannotate <- dat.meta %>%
    left_join(., dat.LL.best)
  
  m1 <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = ctype.from.LL)) +
    geom_point() + 
    theme_bw() + 
    facet_wrap(~ctype.from.LL) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1)
  
  m2 <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2)
  
  # Write outputs -----------------------------------------------------------
  
  fwrite(dat.meta.reannotate, file = outmeta, sep = "\t")
  dev.off()
  
  
}



# 
# # annotate some of the celltypes
# 
# dat.meta.to.change <- subset(dat.meta.reannotate, ctype %in% c("AllCells", "IL7RLinNeg", "LinNeg", "LSK", "HSPCs"))  %>%
#   dplyr::mutate(ctype = ctype.from.LL)
# dat.meta.keep <- subset(dat.meta.reannotate, !ctype %in% c("AllCells", "IL7RLinNeg", "LinNeg", "LSK", "HSPCs"))
# 
# dat.meta.changed <- rbind(dat.meta.to.change, dat.meta.keep)
# 
# m2 <- ggplot(dat.meta.changed, aes(x = umap1, y = umap2, color = ctype)) +
#   geom_point() + 
#   theme_bw() + 
#   facet_wrap(~ctype) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# JFuncs::multiplot(m1, m2, cols = 2)







