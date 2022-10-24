# Jake Yeung
# Date of Creation: 2022-07-24
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/4-get_dynamic_bins_by_mark.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k4me1", "k9me3", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks



# Load count mats ---------------------------------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/processed_data"

mats.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(indir, paste0("countmat_bonemarrow_allmerged_", jmark, ".txt.gz"))
  mat <- data.frame(fread(inf))
  rownames(mat) <- mat$V1
  mat$V1 <- NULL
  return(mat)
})


# Get bins ----------------------------------------------------------------

bins.lst <- lapply(mats.lst, function(jmat){
  rownames(jmat)
})

print(bins.lst[[1]][[1]])

# Write output ------------------------------------------------------------

dat.bins.lst <- lapply(bins.lst, function(bvec){
  chromos <- sapply(bvec, function(b) strsplit(b, ":")[[1]][[1]])
  chromosnoprefix <- gsub("^chr", "", chromos)
  startends <- sapply(bvec, function(b) strsplit(b, ":")[[1]][[2]]) 
  starts <- sapply(startends, function(se) strsplit(se, "-")[[1]][[1]])
  ends <- sapply(startends, function(se) strsplit(se, "-")[[1]][[2]])
  dat <- data.frame(Chr = chromosnoprefix, Start = starts, Ends = ends, bname = bvec, stringsAsFactors = FALSE)
})

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/dynamic_bins"
for (jmark in jmarks){
  print(jmark)
  outf.tmp <- file.path(outdir, paste0("dynamic_bins_50kb.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(dat.bins.lst[[jmark]], file = outf.tmp, sep = "\t", col.names = FALSE)
}


