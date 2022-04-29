# Jake Yeung
# Date of Creation: 2022-04-18
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/20-clean_up_metadata.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k4me3", "k4me1", "k9me3", "k27me3"); names(jmarks) <- jmarks

# Load metas --------------------------------------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(indir, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt"))
  fread(inf)
})

dat.meta.colors <- subset(dat.meta.lst[[1]], select = c(ctype.from.LL, colcode))
dat.meta.colors <- dat.meta.colors[duplicated(dat.meta.colors), ]


# Assign Tcells  ----------------------------------------------------------

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

ggplot(dat.meta.lst$k4me1, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~ctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

