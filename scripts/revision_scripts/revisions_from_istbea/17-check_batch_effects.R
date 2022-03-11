# Jake Yeung
# Date of Creation: 2022-03-01
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/17-check_batch_effects.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_repressive_cleaned/metadata_cleaned.k27me3.2022-02-16.txt"
# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_repressive_cleaned/metadata_cleaned.k27me3.2022-02-16.txt"
# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_repressive_cleaned_Tcells/metadata_reannotate_from_LLmat_dynamicbins.k27me3.2022-02-23.txt"
# inf.meta <- 
dat.meta <- fread(inf.meta) 

dat.meta <- dat.meta %>%
  rowwise() %>%
  mutate(is.pdc.ll = ctype.from.LL == "pDCs", 
         is.pdc = ctype == "pDCs")

ggplot(dat.meta, aes(x = umap1, y = umap2, color = is.pdc)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~ctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_identity() + 
  facet_wrap(~ctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta, aes(x = umap1, y = umap2, color = is.pdc.ll)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
