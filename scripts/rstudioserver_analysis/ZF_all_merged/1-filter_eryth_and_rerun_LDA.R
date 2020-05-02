# Jake Yeung
# Date of Creation: 2020-04-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/1-filter_eryth_and_rerun_LDA.R
# Remove eryth and rerun LDA 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

# Load LDA ----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000")

count.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  infname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
  inf <- file.path(inmain, infname)
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  return(count.mat)
})



# Load clustering table ---------------------------------------------------

dat.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.annot <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  dat.annot <- fread(inf.annot)
  dat.annot$mark <- jmark
  return(dat.annot)
})

# Remove eryth  -----------------------------------------------------------

eryth.cells.lst <- lapply(jmarks, function(jmark){
  eryth.cells <- subset(dat.annot.lst[[jmark]], grepl("eryth", cluster))$cell
  return(eryth.cells)
})


# Write table -------------------------------------------------------------

count.mat.filt.lst <- lapply(jmarks, function(jmark){
  cells.remove <- eryth.cells.lst[[jmark]]
  cells.all <- colnames(count.mat.lst[[jmark]])
  cells.keep <- !colnames(count.mat.lst[[jmark]]) %in% cells.remove
  count.mat.filt <- count.mat.lst[[jmark]][, cells.keep]
  return(count.mat.filt)
})

lapply(count.mat.lst, dim)
lapply(count.mat.filt.lst, dim)

# Write new tables to output ----------------------------------------------

outdir <- file.path(hubprefix, "/jyeung/data/zebrafish_scchic/from_rstudio/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.corrfilt.eryth_removed")

for (jmark in jmarks){
  print(jmark)
  saveRDS(count.mat.filt.lst[[jmark]], file = file.path(outdir, paste0("count_mat.", jmark, ".eryth_filt.", Sys.Date(), ".rds")))
}

