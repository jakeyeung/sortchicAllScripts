# Jake Yeung
# Date of Creation: 2020-02-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/6b-split_KeepMorePlates_to_unenriched.R
# KeepMorePlates is missing "Unenriched" matrix

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM.final_Unenriched_and_AllMerged"
dir.create(outdir)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

whitelst <- c("B6-13W1-BM-H3", "PZ-ChIC-Bl6-H3", "PZ-ChIC-Bl6-BM-H3K4me1-Index")
whitelst.grep <- paste(whitelst, collapse = "|")

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
plates.whitelist <- lapply(plates.lst, function(p){
  p[grepl(whitelst.grep, p)]
}) %>%
  unlist()



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
  cells.keep.i <- which(plates.all %in% plates.whitelist)
  return(list(AllMerged = count.mat, Unenriched = count.mat[, cells.keep.i]))
})

# check plates

plates.all <- lapply(count.mats.all, function(x){
  plates.merged <- unique(sapply(colnames(x$AllMerged), function(x) ClipLast(x, jsep = "_")))
  plates.unenriched <- unique(sapply(colnames(x$Unenriched), function(x) ClipLast(x, jsep = "_")))
  return(list(AllMerged = plates.merged, Unenriched = plates.unenriched))
})

lapply(plates.all, function(x) print(x$Unenriched))
lapply(plates.all, function(x) print(x$AllMerged))


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



