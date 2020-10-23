# Jake Yeung
# Date of Creation: 2020-08-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/K562_analysis/merge_old_and_new_for_LDA.R
# Merge old and new 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load  -------------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cleaned_count_tables_for_lda_and_projections.top5000.intersect"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cleaned_count_tables_for_lda_and_projections"
outdir <- file.path(indir, "merged_tables")
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  inf.old <- file.path(indir, paste0("count_mats_old_binsize_50000_genomewide.", jmark, ".old.rds"))
  inf.new <- file.path(indir, paste0("count_mats_old_binsize_50000_genomewide.", jmark, ".new.rds"))
  outf <- file.path(outdir, paste0("count_mats_old_binsize_50000_genomewide.", jmark, ".merged.rds"))
  if (file.exists(outf)){
    print(paste("Skipping:", outf))
    next
  }
  mat.old <- readRDS(inf.old)
  mat.new <- readRDS(inf.new)
  assertthat::assert_that(identical(rownames(mat.old), rownames(mat.new)))
  mat.merged <- cbind(mat.old, mat.new)
  # save output
  saveRDS(object = mat.merged, file = outf)
}
