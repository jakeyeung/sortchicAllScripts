# Jake Yeung
# Date of Creation: 2020-09-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/filt_cells_after_LDA_by_var.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

# Load LDA  ---------------------------------------------------------------

var.cutoff <- 1.5
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

hubprefix <- "/home/jyeung/hub_oudenaarden"
inmain <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046")
jsuffix <- "_counts_filt.2020-09-12.K-30"
tm.result.lst <- lapply(jmarks, function(jmark){
  dname <- paste0("lda_outputs.count_mat_", jmark, jsuffix, ".binarize.FALSE")
  indir <- file.path(inmain, dname)
  fname <- paste0("ldaOut.count_mat_", jmark, "_counts_filt.2020-09-12.K-30.Robj")
  inf.lda <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.lda))
  print(inf.lda)
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(tm.result)
})

dat.var.lst <- lapply(jmarks, function(jmark){
  tm.result <- tm.result.lst[[jmark]]
  dat.impute <- t(log2(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute, jchromos)
  dat.var$mark <- jmark
  return(dat.var)
}) 


dat.var.long <- bind_rows(dat.var.lst)

# Get low var cells -------------------------------------------------------

bad.cells <- subset(dat.var.long, cell.var.within.sum.norm <= var.cutoff)$cell


# Write new tables?  ------------------------------------------------------


indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.VAN5046"
outdir <- file.path(indir, paste0("varfilt"))
dir.create(outdir)

infs.rds <- list.files(indir, pattern = "*.rds", full.names = TRUE)

for (inf in infs.rds){
  bname <- basename(inf)
  bname <- gsub(pattern = ".rds", replacement = "", x = bname)
  outname <- paste0(bname, ".varfilt_", var.cutoff, ".rds")
  outf <- file.path(outdir, outname)
  mat <- readRDS(inf)
  cells.keep <- colnames(mat)[!colnames(mat) %in% bad.cells]
  
  mat.filt <- mat[, cells.keep]
  
  print("Dim before filt") 
  print(dim(mat))
  print("Dim after filt") 
  print(dim(mat.filt))
  
  saveRDS(mat.filt, file = outf)
}


