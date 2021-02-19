# Jake Yeung
# Date of Creation: 2021-02-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/2-make_dendrograms_from_LDA_pseudobulk.dynamic_bins.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load metadata -----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname))
})

# Load LDA  ---------------------------------------------------------------


jmark <- "H3K4me3"
# jsuffixmain <- paste0("dynamic_bins.50kb")
jsuffixmain <- paste0("dynamic_bins.50kb")
# jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"
# jsuffix <- "celltype_specific_genes.TSS_10kb.txt"
# jsuffix <- "celltype_specific_genes.TSS_TES.txt"
# dat.impute.pbulk.lst <- lapply(jmarks, function(jmark){


out.lst <- lapply(jmarks, function(jmark){
  # jsuffixmain <- "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.H3K9me3.2021-02-15.txt"
  
  if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
    jsuffix <- paste0(jsuffixmain, ".", jmark, ".2021-02-15.txt")
  } else {
    jsuffix <- paste0("dynamic_bins.50kb.", jmark, ".txt")
  }
  
  # jsuffix <- paste0(jsuffixmain, ".", jmark, ".txt")
  print(jmark)
  
  if (jsuffixmain == "DE_bins_all_marks_top_6085_dists_to_TSS.annot_table"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes.top_6085/lda_outputs.count_tables_merged.", jmark, ".", jsuffix, ".K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".", jsuffix, ".K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM.varfilt.dynamic_bins_genes/lda_outputs.count_tables_merged.", jmark, ".", jsuffix, ".K-30.binarize.FALSE/ldaOut.count_tables_merged.", jmark, ".", jsuffix, ".K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda))
  
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.impute <- t(tm.result$topics %*% tm.result$terms)
  
  cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
  dat.impute.pbulk <- do.call(cbind, SumAcrossClusters(dat.impute, cnames.keep.lst))
  dat.impute.pbulk <- sweep(dat.impute.pbulk, MARGIN = 2, STATS = colSums(dat.impute.pbulk), FUN = "/", check.margin = TRUE)
  
  dat.raw.pbulk <- SumAcrossClusters(count.mat, cnames.keep.lst)
  dat.raw.pbulk <- do.call(cbind, dat.raw.pbulk)
  # normalize
  dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * 1000 + 1
  return(list(dat.impute.pbulk = dat.impute.pbulk, dat.raw.pbulk = dat.raw.pbulk))
})

dat.impute.pbulk.lst <- lapply(out.lst, function(x) x$dat.impute.pbulk)
dat.raw.pbulk.lst <- lapply(out.lst, function(x) x$dat.raw.pbulk)

# jmethod <- "complete"
jmethod <- "ward.D2"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"
outpdf <- file.path(outdir, paste0("heatmaps.", jsuffixmain, ".", jmethod, ".pdf"))

pdf(outpdf, useDingbats = FALSE)
for (jmark in jmarks){
  print(jmark)
  heatmap3::heatmap3(log2(dat.impute.pbulk.lst[[jmark]]), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste(jmark, jmethod, jsuffixmain), margins = c(5, 8), cexRow = 0.5, method = jmethod)
}
for (jmark in jmarks){
  print(jmark)
  outtxt <- file.path(outdir, paste0("table.", jsuffixmain, ".", jmethod, ".txt"))
  write.table(x = log2(dat.raw.pbulk.lst[[jmark]]), file = outtxt, row.names = TRUE, col.names = NA)
  heatmap3::heatmap3(log2(dat.raw.pbulk.lst[[jmark]]), Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, main = paste("raw", jmark, jmethod, jsuffixmain), margins = c(5, 8), cexRow = 0.5, method = jmethod)
}
dev.off()



