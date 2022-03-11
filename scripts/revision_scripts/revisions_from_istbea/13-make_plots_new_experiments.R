# Jake Yeung
# Date of Creation: 2022-02-11
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/13-make_plots_new_experiments.R
# 

rm(list=ls()) 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

# Load meta ---------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/rename"
infs.meta <- lapply(jmarks, function(jmark){
  file.path(indir.meta, paste0("metadata_", jmark, ".txt"))
})

dat.meta.lst <- lapply(infs.meta, function(jinf){
  fread(jinf) %>%
    filter(ctype.from.LL != "Tcells")
})

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    facet_wrap(~ctype.from.LL) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
})

m.louvain.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = ctype.from.LL)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    facet_wrap(~louvain) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
})

JFuncs::multiplot(m.lst$k9me3, m.louvain.lst$k9me3, cols = 2)
JFuncs::multiplot(m.lst$k4me3, m.louvain.lst$k4me3, cols = 2)
JFuncs::multiplot(m.lst$k4me1, m.louvain.lst$k4me1, cols = 2)
JFuncs::multiplot(m.lst$k27me3, m.louvain.lst$k27me3, cols = 2)

# Assign H3K9me3 using louvain  -------------------------------------------


dat.k9.new <- dat.meta.lst$k9me3 %>%
  
  

# Plot final UMAPs --------------------------------------------------------




# Show specific cell types (old and groundtruth) --------------------------





