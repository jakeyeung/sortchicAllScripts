# Jake Yeung
# Date of Creation: 2022-07-26
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/1-add_TSS_to_k9me3_cluster_specific_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf.bins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/H3K9me3_H3K4me1_analysis/H3K9me3_cluster_specific_bins.keeptop_500.2022-07-22.txt"
inf.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/refseq/MmRefseqTss.chromorenamed.10000.again.nochromo.sorted.bed"

dat.bins <- fread(inf.bins, header = FALSE)
dat.tss <- fread(inf.tss, header = FALSE)
dat.tss$V5 <- NULL

# add coord to inf.tss
tss.coords <- paste(dat.tss$V1, paste(dat.tss$V2, dat.tss$V3, sep = "-"), sep = ":")
tss.orig <- dat.tss$V4
dat.tss$V4 <- paste(tss.coords, tss.orig, sep = ";")

# Combine -----------------------------------------------------------------

dat.merged <- rbind(dat.bins, dat.tss)


# Save output -------------------------------------------------------------

fout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/refseq/MmRefseqTss.chromorenamed.10000.with_k9_bins_keeptop_500.", Sys.Date(), ".txt")

fwrite(dat.merged, file = fout, sep = "\t", quote = FALSE, col.names = FALSE)


