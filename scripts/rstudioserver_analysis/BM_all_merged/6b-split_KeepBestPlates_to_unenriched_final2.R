# Jake Yeung
# Date of Creation: 2020-02-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/6b-split_KeepBestPlates_to_unenriched_final2.R
# Keep best plates not just keep more plates

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM.KeepBestPlates.final_Unenriched_and_AllMerged"
dir.create(outdir)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

whitelst <- c("B6-13W1-BM-H3", "PZ-ChIC-Bl6-H3", "PZ-ChIC-Bl6-BM-H3K4me1-Index")
whitelst.grep <- paste(whitelst, collapse = "|")
blacklst <- c("B6-13W1-BM-H3K4me3-4", "B6-13W1-BM-H3K27me3-4", "B6-13W1-BM-H3K9me3-4")
blacklst.grep <- paste(blacklst, collapse = "|")

whitelst.sorted <- c("-stem-cells-", "-lin-", "Linneg", "BMSC")
whitelst.sorted.grep <- paste(whitelst.sorted, collapse = "|")


# for (jmark in jmarks){
plates.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  # H3K4me1 is special: Linneg added later. Antibody changed halfway through. There are bad *plates* in the new whole bonemarrow.
  if (jmark == "H3K4me1"){
    inf <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM_includeNewLinnegH3K4me1.2020-01-31/BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.WithLinneg.rds")
  } else {
    inf <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM.keepPlates/BM_", jmark, ".varcutoff_0.3.KeepAllPlates.rds")
  } 
  print(paste("Loading inf:", inf))
  count.mat <- readRDS(inf)
  plates <- unique(sapply(colnames(count.mat), function(x) ClipLast(x, jsep = "_")))
  return(plates)
})

# check
plates.whitelist.unenriched <- lapply(plates.lst, function(p){
  p[grepl(whitelst.grep, p) & !grepl(blacklst.grep, p)]
}) %>%
  unlist()

plates.whitelist.sorted <- lapply(plates.lst, function(p){
  p[grepl(whitelst.sorted.grep, p)]
}) %>%
  unlist()

plates.whitelist.allmerged <- c(plates.whitelist.unenriched, plates.whitelist.sorted)


# Write AllMerged and Unenriched ------------------------------------------

count.mats.all <- lapply(jmarks, function(jmark){
  print(jmark)
  # H3K4me1 is special: Linneg added later. Antibody changed halfway through. There are bad *plates* in the new whole bonemarrow.
  if (jmark == "H3K4me1"){
    inf <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM_includeNewLinnegH3K4me1.2020-01-31/BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.WithLinneg.rds")
  } else {
    inf <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM.keepPlates/BM_", jmark, ".varcutoff_0.3.KeepAllPlates.rds")
  } 
  print(paste("Loading inf:", inf))
  count.mat <- readRDS(inf)
  plates.all <- sapply(colnames(count.mat), function(x) ClipLast(x, jsep = "_"))
  cells.all <- colnames(count.mat)
  cells.keep.i <- which(plates.all %in% plates.whitelist.unenriched)
  cells.keep.i.allmerged <- which(plates.all %in% plates.whitelist.allmerged)
  return(list(AllMerged = count.mat[, cells.keep.i.allmerged], Unenriched = count.mat[, cells.keep.i]))
})

# check plates

plates.all <- lapply(count.mats.all, function(x){
  plates.merged <- unique(sapply(colnames(x$AllMerged), function(x) ClipLast(x, jsep = "_")))
  plates.unenriched <- unique(sapply(colnames(x$Unenriched), function(x) ClipLast(x, jsep = "_")))
  return(list(AllMerged = plates.merged, Unenriched = plates.unenriched))
})

lapply(plates.all, function(x) x$Unenriched)
lapply(plates.all, function(x) x$AllMerged)


# Write count mat to output  ----------------------------------------------

lapply(jmarks, function(jmark){
  print(jmark)
  mat.lst <- count.mats.all[[jmark]]
  outname <- paste0("BM_", jmark, "_varfilt_countmat.", Sys.Date())
  outprefix <- file.path(outdir, outname)
  count.mat.allmerged <- mat.lst$AllMerged
  count.mat.unenriched <- mat.lst$Unenriched
  outf.AllMerged <- paste0(outprefix, ".AllMerged.rds")
  outf.Unenriched <- paste0(outprefix, ".Unenriched.rds")
  
  print(paste("Saving AllMerged of", paste(dim(count.mat.allmerged), collapse = ","), "to", outf.AllMerged))
  saveRDS(mat.lst$AllMerged, file = outf.AllMerged)
  print(paste("Saving Unenriched of", paste(dim(count.mat.unenriched), collapse = ","), "to", outf.Unenriched))
  saveRDS(mat.lst$Unenriched, file = outf.Unenriched)
})



