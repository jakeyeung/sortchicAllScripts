# Jake Yeung
# Date of Creation: 2020-11-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/1-QC_plots_K562.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(ggrastr)

library(topicmodels)

library(JFuncs)
library(scchicFuncs)

# Load K562 for QC plots --------------------------------------------------

jchromos <- paste("chr", c(seq(22), "X", "Y"), sep = "")
jchromos.nochr <- paste(c(seq(22), "X", "Y"), sep = "")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jdate <- "2020-11-30"
indir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2")
intxt <- file.path(indir, paste0("K562_QC_plots.", jdate, ".add_fracnonzeros.txt"))

dats.raw <- fread(intxt)

dats <- lapply(jmarks, function(jmark){
  print(jmark)
  intxt.tmp <- file.path(indir, paste0("K562_QC_plots.2020-11-30.add_fracnonzeros.", jmark, ".good_cells.txt"))
  fread(intxt.tmp)
})

dats.long <- dats %>%
  bind_rows()

