# Jake Yeung
# Date of Creation: 2020-11-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/5-analyze_celltype_specific_peaks.R
# Reproduce celltype specific peaks? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load LDA from peaks  ----------------------------------------------------


