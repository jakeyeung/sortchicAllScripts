# Jake Yeung
# Date of Creation: 2020-04-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/test_collapse_gene_by_row_TSS.R
# Test that selecting best TSS is working as expected

rm(list=ls())

library(dplyr)
library(Matrix)
library(scchicFuncs)

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS/H3K4me1.countTableTSS.mapq_40.TSS_10000.blfiltered.rds"

mat <- readRDS(inf)

mat[1:5, 1:5]

mat.collapsed <- CollapseRowsByGene(mat, as.long=FALSE)

# check Sox17

mat.filt <- mat[grepl("Sox17", rownames(mat)), ]

mat.filt.collapsed <- mat.collapsed[grepl("Sox17", rownames(mat.collapsed)), ]

sum(mat.filt.collapsed)
rowSums(mat.filt)


