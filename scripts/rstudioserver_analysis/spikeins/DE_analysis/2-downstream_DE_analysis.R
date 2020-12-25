# Jake Yeung
# Date of Creation: 2020-12-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/2-downstream_DE_analysis.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs) 



# Check DEs ---------------------------------------------------------------


inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_bins.H3K4me1.2020-12-12.newannot2.witherrors.MoreBins.RData"
load(inf, v=T)

# get distribution of log2FCs 




