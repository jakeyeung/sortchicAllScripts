# Jake Yeung
# Date of Creation: 2020-11-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/3-load_count_tables_fill_zeros_for_LDA.TSS.10kb.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)
library(scchicFuncs)


jsuffix <- "from_TSS"
jdist <- "10000"
hubprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.", jsuffix, ".dist_", jdist))
dir.create(outdir)

# indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_", jsuffix))
indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_TSS/dist_10000"))
assertthat::assert_that(dir.exists(indir))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

mats.bymark.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  infs.bymark <- list.files(path = indir, pattern = paste0("BM_round1_round2_merged_", jmark, ".*.txt"), full.names = TRUE)
  names(infs.bymark) <- sapply(infs.bymark, basename)
  
  mats.bymark <- lapply(infs.bymark, function(inf){
    ReadMatTSSFormat(inf, add.coord = TRUE, sort.rnames = TRUE)
  })
  
  jall.rnames <- sort(unique(unlist(lapply(mats.bymark, rownames))))
  mats.bymark.fill <- cbind.fill.lst(mats.bymark, all.rnames = jall.rnames, fill = 0)
  return(mats.bymark.fill)
})


# Write to output ---------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_", jsuffix, ".", jmark, ".dist_", jdist, ".rds"))
  assertthat::assert_that(!file.exists(outf))
  saveRDS(mats.bymark.lst[[jmark]], file = outf)
}




