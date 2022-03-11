# Jake Yeung
# Date of Creation: 2022-01-26
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/1-get_dynamic_bins_four_marks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Get dynamic bins of the four marks --------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarks

markold2new <- hash::hash(jmarks, jmarksnew)

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/count_tables.BM.dynamic_bins_TSS_TES_regions"
outmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed_dynamic_bins"

infs.meta <- lapply(jmarks, function(jmark){
  file.path(indir, paste0("dynamic_bins.50kb.corrected_DE_tables.", jmark, ".2021-04-07.txt"))
})

bins.keep.lst <- lapply(infs.meta, function(jinf){
  fread(jinf, col.names = c("chromo", "start", "end", "rname"))
})


# Load mats ---------------------------------------------------------------

inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed"

infs.lst.lst <- lapply(jmarksnew, function(jmark){
  dname <- paste0("BM_", jmark)
  indir <- file.path(inmain, dname)
  infs.mat <- list.files(path = indir, pattern = "count_mat.*.rds", full.names = TRUE)
  names(infs.mat) <- sapply(infs.mat, function(jinf) strsplit(basename(jinf), split = "\\.")[[1]][[1]])
  return(infs.mat)
})


mats.lst.lst <- lapply(infs.lst.lst, function(jinf.lst){
  lapply(jinf.lst, function(jinf){
    readRDS(jinf)
  })
})


# Use it to filter counts  ------------------------------------------------

mats.filt.lst.lst <- lapply(jmarks, function(jmark){
  jmat.lst <- mats.lst.lst[[jmark]]
  bins.keep <- bins.keep.lst[[jmark]]$rname
  lapply(jmat.lst, function(jmat){
    jmat.filt <- jmat[bins.keep, ]
    assertthat::assert_that(nrow(jmat.filt) > 0)
    return(jmat.filt)
  })
})


# Check empty rows --------------------------------------------------------



# Write new counts --------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  jmarknew.tmp <- markold2new[[jmark]]
  outdir <- file.path(outmain, paste0("BM_", jmarknew.tmp))
  dir.create(outdir)
  jmat.lst <- mats.filt.lst.lst[[jmark]]
  jnames <- names(jmat.lst)
  for (jname in jnames){
    print(jname)
    outf.tmp <- file.path(outdir, paste0(jname, ".", jmarknew.tmp, ".", Sys.Date(), ".rds"))
    mat.tmp <- mats.filt.lst.lst[[jmark]][[jname]]
    print("Dim before")
    print(dim(mat.tmp))
    
    # check duplicated cells
    cells <- colnames(mat.tmp)
    if (length(cells) != length(unique(cells))){
      print("Removing duplicated cells")
      cells.duplicated <- duplicated(colnames(mat.tmp))
      mat.tmp <- mat.tmp[, !cells.duplicated]
      print(dim(mat.tmp))
    }
    
    bins.empty <- rowSums(mat.tmp) == 0
    cells.empty <- colSums(mat.tmp) == 0
    print(paste("Empty bins:", length(which(bins.empty))))
    print(paste("Empty cells:", length(which(cells.empty))))
    if (length(which(bins.empty)) > 0){
      print("Found empty bins, removing...")
      mat.tmp <- mat.tmp[!bins.empty, ]
      print("Dim after removing empty bins")
      print(dim(mat.tmp))
    }
    if (length(which(cells.empty)) > 0){
      print("Found empty cells, removing...")
      mat.tmp <- mat.tmp[, !cells.empty]
      print("Dim after removing empty cells")
      print(dim(mat.tmp))
    }
    saveRDS(mat.tmp, file = outf.tmp)
  }
}

