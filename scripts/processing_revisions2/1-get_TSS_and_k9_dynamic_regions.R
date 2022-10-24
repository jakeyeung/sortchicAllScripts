# Jake Yeung
# Date of Creation: 2022-07-18
# File: ~/projects/scchic/scripts/revisions2/1-get_TSS_and_k9_dynamic_regions.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


bsize <- 10000
infbed <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/refseq/MmRefseqTss.chromorenamed.", bsize, ".again.nochromo.sorted.bed")
datbed <- fread(infbed, header = FALSE)

print(head(datbed))

# Load K9me3 bins  -------------------------------------------------------------

inf.k9 <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/processed_data/countmat_bonemarrow_allmerged_k9me3.txt.gz"
mat.k9 <- fread(inf.k9)
bins.keep <- mat.k9$V1
bins.keep.nochromo <- sapply(bins.keep, function(b){
  gsub("^chr", "", b)
})
chromos <- sapply(bins.keep.nochromo, function(b) strsplit(b, ":")[[1]][[1]])
startends <- sapply(bins.keep.nochromo, function(b) strsplit(b, ":")[[1]][[2]])
starts <- sapply(startends, function(x) strsplit(x, "-")[[1]][[1]])
ends <- sapply(startends, function(x) strsplit(x, "-")[[1]][[2]])

datbed.k9 <- data.frame(V1 = chromos, V2 = starts, V3 = ends, V4 = bins.keep.nochromo, V5 = "+", stringsAsFactors = FALSE)
rownames(datbed.k9) <- NULL

head(datbed.k9)

# Write new bed  ----------------------------------------------------------

datbednew <- rbind(datbed, datbed.k9)


# write output ------------------------------------------------------------

outf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/refseq/MmRefseqTss.chromorenamed.", bsize, ".with_k9_bins.", Sys.Date(), ".txt")
fwrite(datbednew, file = outf, quote = FALSE, sep = "\t", col.names = FALSE)
