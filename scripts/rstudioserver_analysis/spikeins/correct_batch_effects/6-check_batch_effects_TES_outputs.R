# Jake Yeung
# Date of Creation: 2021-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/6-check_batch_effects_TES_outputs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)



hubprefix <- "/home/jyeung/hub_oudenaarden"

inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/robjs_batch_correction_output/batch_corrected_imputed_values.TES.k4_k9.mat.namesfix.2021-02-12.TES.from_LDA.RData"))
load(inf, v=T)

