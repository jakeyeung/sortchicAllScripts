# Jake Yeung
# Date of Creation: 2020-08-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/quality_control_K562_spikein_recovery.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Load dat  ---------------------------------------------------------------


# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters.filt_spikeins"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters"

infs.chromo <- list.files(indir, pattern = "*.csv", all.files = TRUE, full.names = TRUE)
jnames <- sapply(infs.chromo, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]], simplify = TRUE, USE.NAMES = FALSE)
names(jnames) <- jnames
names(infs.chromo) <- jnames

jchromos <- paste(c(seq(19), "X", "Y"), sep = "")

jdats <- lapply(jnames, function(jname){
  jinf <- infs.chromo[[jname]]
  jdat <- scchicFuncs::GetChromoCounts(inf = jinf, spikeinchromo = "J02459.1", chromos.keep = jchromos)
  jdat$experi <- jname
  jdat <- subset(jdat, chromo == 1)
  return(jdat)
}) %>%
  bind_rows()

ggplot(jdats, aes(x = experi, y = log2(spikeincounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  







