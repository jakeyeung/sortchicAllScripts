# Jake Yeung
# Date of Creation: 2020-04-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/check_DE_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(Seurat)

library(JFuncs)
library(scchicFuncs)



# Load DE  ----------------------------------------------------------------

inf <- "/home/jyeung/data/from_rstudioserver/zebrafish.2020-04-26/diff_exprs_Chloe_seurat.full.rds"
de.out <- readRDS(inf)

plot(density(log10(abs(de.out$avg_logFC) + 1)))

