# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/get_mean_exprs_across_bins.R
# Mean exprs

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(JFuncs)
library(scchicFuncs)


# inf1 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf2 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf3 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf4 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K27me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K27me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf5 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf6 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf7 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf8 <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K27me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K27me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jprefixs <- c("Unenriched", "AllMerged")

inf.gc <- "/home/jyeung/hpc/scChiC/from_rstudioserver/rdata_robjs/gr_gc_dat.RData"
load(inf.gc, v=T)

outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/gc_signal_across_plates"
assertthat::assert_that(dir.exists(outdir))


for (jprefix in jprefixs){
  print(jprefix)
  for (jmark in jmarks){
    print(jmark)
    inf <- paste0("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_", jprefix, "_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_", jprefix, "_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj")
    assertthat::assert_that(file.exists(inf))
    load(inf, v=T)
    
    count.mat.long <- data.table::melt(as.matrix(count.mat))
    colnames(count.mat.long) <- c("coord", "cell", "count")
    # count.mat.long <- count.mat.long %>%
    #   group_by(cell) %>%
    #   mutate(frac.umi = count / sum(count))
    
    experi.hash <- hash(colnames(count.mat), sapply(colnames(count.mat), ClipLast, jsep = "_"))
    
    experis <- unique(hash::values(experi.hash))
    names(experis) <- experis
    count.mat.sum <- lapply(experis, function(experi){
      print(experi)
      cells <- grep(experi, colnames(count.mat), value = TRUE)
      count.mat.sum <- data.frame(coord = rownames(count.mat), experi = experi, count = rowSums(count.mat[, cells]))
      return(count.mat.sum)
    }) %>%
      bind_rows()
    
    
    # count.mat.long$experi <- sapply(as.character(count.mat.long$cell), function(x) ClipLast(x, jsep = "_"))
    # count.mat.long$experi <- sapply(as.character(count.mat.long$cell), function(x) ClipLast(x, jsep = "_"))
    count.mat.long$experi <- sapply(as.character(count.mat.long$cell), function(x) experi.hash[[x]])
    
    
    # count.mat.sum <- count.mat.long %>%
    #   group_by(coord, experi) %>%
    #   summarise(count = mean(count))
    
    gr.gc.dat$start <- as.integer(gr.gc.dat$start)
    gr.gc.dat$end <- as.integer(gr.gc.dat$end)
    gr.gc.dat <- gr.gc.dat %>%
      rowwise() %>%
      mutate(coord = paste(chromo, paste(start, end, sep = "-"), sep = ":"))
    
    gr.gc.dat.merge <- left_join(gr.gc.dat, count.mat.sum)
    
    gr.gc.dat.merge.subset <- gr.gc.dat.merge %>%
      group_by(experi) %>%
      sample_frac(tbl = ., size = 0.15, replace = FALSE)
    
    pdf(file.path(outdir, paste0("gc_count_plots.", jprefix, ".", jmark, ".pdf")), useDingbats = FALSE)
      ggplot(gr.gc.dat.merge.subset, aes(x = gc, y = count)) + geom_point(alpha = 0.1) + facet_wrap(~experi, ncol = 1) + geom_smooth(method = "lm")
      ggplot(gr.gc.dat.merge.subset, aes(x = gc, y = count)) + geom_point(alpha = 0.1) + facet_wrap(~experi) + geom_smooth(method = "lm")
      ggplot(gr.gc.dat.merge.subset, aes(x = gc, y = count)) + geom_point(alpha = 0.1) + facet_wrap(~experi) + geom_smooth(method = "loess")
    dev.off()
  }
}


# for (inf in infs.lst){
#   assertthat::assert_that(file.exists(inf))
#   load(inf, v=T)
#   
#   inf.gc <- "/home/jyeung/hpc/scChiC/from_rstudioserver/rdata_robjs/gr_gc_dat.RData"
#   load(inf.gc, v=T)
#   
#   
#   count.mat.long <- data.table::melt(as.matrix(count.mat))
#   colnames(count.mat.long) <- c("coord", "cell", "count")
#   # count.mat.long <- count.mat.long %>%
#   #   group_by(cell) %>%
#   #   mutate(frac.umi = count / sum(count))
#   
#   experi.hash <- hash(colnames(count.mat), sapply(colnames(count.mat), ClipLast, jsep = "_"))
#   
#   experis <- unique(hash::values(experi.hash))
#   names(experis) <- experis
#   count.mat.sum <- lapply(experis, function(experi){
#     print(experi)
#     cells <- grep(experi, colnames(count.mat), value = TRUE)
#     count.mat.sum <- data.frame(coord = rownames(count.mat), experi = experi, count = rowSums(count.mat[, cells]))
#     return(count.mat.sum)
#   }) %>%
#     bind_rows()
#   
#   
#   # count.mat.long$experi <- sapply(as.character(count.mat.long$cell), function(x) ClipLast(x, jsep = "_"))
#   # count.mat.long$experi <- sapply(as.character(count.mat.long$cell), function(x) ClipLast(x, jsep = "_"))
#   count.mat.long$experi <- sapply(as.character(count.mat.long$cell), function(x) experi.hash[[x]])
#   
#   
#   # count.mat.sum <- count.mat.long %>%
#   #   group_by(coord, experi) %>%
#   #   summarise(count = mean(count))
#   
#   gr.gc.dat$start <- as.integer(gr.gc.dat$start)
#   gr.gc.dat$end <- as.integer(gr.gc.dat$end)
#   gr.gc.dat <- gr.gc.dat %>%
#     rowwise() %>%
#     mutate(coord = paste(chromo, paste(start, end, sep = "-"), sep = ":"))
#   
#   gr.gc.dat.merge <- left_join(gr.gc.dat, count.mat.sum)
#   
#   gr.gc.dat.merge.subset <- gr.gc.dat.merge %>%
#     group_by(experi) %>%
#     sample_frac(tbl = ., size = 0.15, replace = FALSE)
#   
#   ggplot(gr.gc.dat.merge.subset, aes(x = gc, y = count)) + geom_point(alpha = 0.1) + facet_wrap(~experi, ncol = 1) + geom_smooth(method = "lm")
#   ggplot(gr.gc.dat.merge.subset, aes(x = gc, y = count)) + geom_point(alpha = 0.1) + facet_wrap(~experi) + geom_smooth(method = "lm")
#   ggplot(gr.gc.dat.merge.subset, aes(x = gc, y = count)) + geom_point(alpha = 0.1) + facet_wrap(~experi) + geom_smooth(method = "loess")
#   
# }
# 
# 
# 