# Jake Yeung
# Date of Creation: 2021-02-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/19-make_heatmap_K4me1_and_K9me3_dynamic_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

reconvert <- TRUE
jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks

# Load batch corrected mat  -----------------------------------------------

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2021-02-05.H3K4me1_H3K9me3.TES_50kb_together.RData"
load(inf, v=T)

if (reconvert){
  dat.adj.dedup <- dat.adj.lst2$H3K4me1[!duplicated(dat.adj.lst2$H3K4me1[, c(1,2)]), ]
  mat.adj.dedup <- data.table::dcast(dat.adj.dedup, formula = rname ~ cell, value.var = "log2exprsadj")
  outrds.tmp <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2021-02-05.H3K4me1_H3K9me3.TES_50kb_together.H3K4me1.rds"
  saveRDS(mat.adj.dedup, file = outrdstmp)
}
 
# 
# # Convert dat adj lst to wide  --------------------------------------------
# 
# mat.adj.lst$H3K4me1 <- mat.adj.dedup
# 
# 
# # Get dat meats -----------------------------------------------------------
# 
# dat.metas <- lapply(jmarks, function(jmark){
#   inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
#   fread(inf.tmp)
# })
# 
# cells.keep.lst <- lapply(jmarks, function(jmark){
#   dat.metas[[jmark]]$cell
# })
# 
# 
# 
# # Order bins?  -------------------------------------------------------------
# 
# 
# 
# 
# 
# # Make heatmap  -----------------------------------------------------------
# 
# jmat.wins.lst <- lapply(jmarks, function(jmark){
#   jmat <- mat.adj.lst[[jmark]][bins.keep.common, cells.keep.lst[[jmark]]]
#   jmat <- apply(jmat, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))
#   jmat <- t(apply(jmat, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
#   return(jmat)
# })
# 
# colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)
# 
# print("Making heatmaps")
# for (jmark in jmarks){
#   heatmap3::heatmap3(jmat.wins.lst[[jmark]], Rowv = TRUE, Colv = NA, scale = "row", ColSideColors = dat.metas[[jmark]]$clustercol,  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8), col = colpalette)
# }
# 
# 
# 
# 
# 
# 
# 
