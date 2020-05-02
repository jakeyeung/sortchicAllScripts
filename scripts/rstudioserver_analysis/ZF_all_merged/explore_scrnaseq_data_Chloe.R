# Jake Yeung
# Date of Creation: 2020-04-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/explore_scrnaseq_data_Chloe.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(Seurat)

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset_PenalizedNegBinRegress.rds")

dat <- readRDS(inf)

# Pseudobulks?  -----------------------------------------------------------




