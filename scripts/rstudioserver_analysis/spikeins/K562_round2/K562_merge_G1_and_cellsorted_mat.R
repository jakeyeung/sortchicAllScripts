# Jake Yeung
# Date of Creation: 2020-09-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/K562_merge_G1_and_cellsorted_mat.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Merge across histone marks  ---------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2"

dat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf1 <- file.path(indir, paste0("K562_count_tables_50000.", jmark, ".cellcyclefilt.rds"))
  inf2 <- file.path(indir, paste0("K562_count_tables_50000.", jmark, ".G1filt.rds"))
  dat1 <- readRDS(inf1)
  dat2 <- readRDS(inf2)
  rnames.keep <- intersect(rownames(dat1), rownames(dat2))
  dat.merged <- cbind(dat1[rnames.keep, ], dat2[rnames.keep, ])
})


# write output
lapply(jmarks, function(jmark){
  print(jmark)
  outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/merged/K562_count_tables_50000.", jmark, ".AllMerged.rds")
  saveRDS(dat.lst[[jmark]], file = outf)
})




