# Jake Yeung
# Date of Creation: 2021-01-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/5-make_heatmaps_K27me3_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)



# Load TSS  ---------------------------------------------------------------





