# Jake Yeung
# Date of Creation: 2019-11-25
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/6-analyze_unenriched_linneg_merged_analysis_oldpipeline.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(Matrix)
library(topicmodels)

library(scchicFuncs)

inf1 <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"
inf2 <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K9me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"
inf3 <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K27me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"
inf4 <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me1_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"

infs <- c(inf1, inf2, inf3, inf4)

pdf("/Users/yeung/data/scchic/pdfs/linneg_analysis.redo/linneg_analysis_all_marks.pdf", useDingbats = FALSE)
for (inf in infs){
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  if (length(out.lda) > 1){
    out.lda <- out.lda[[1]]
  }
  tm.result <- posterior(out.lda)
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  umap.out <- umap(tm.result$topics, config = jsettings)
  
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell, jsep = "-"))
  
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m0 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_manual(values = cbPalette) + 
    facet_wrap(~experi)
  
  
  # do variance
  jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
  jfac <- 10^6
  jpseudo <- 0
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms) * jfac + jpseudo)
  cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)
  dat.umap.long.var <- left_join(dat.umap.long, cells.var.chromo, by = "cell")
  
  m1 <- ggplot(dat.umap.long.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + 
    facet_wrap(~experi)
  m2 <- dat.umap.long.var %>% 
    ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi, ncol = 1)
  m3 <- dat.umap.long.var %>% 
    ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  print(m1)
  print(m2)
  print(m3)
}

dev.off()


