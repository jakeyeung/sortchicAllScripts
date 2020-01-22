# Jake Yeung
# Date of Creation: 2019-12-13
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/explore_K9me3_batch_effects.R
# Why is there separation in K9me3? We didn't see it before? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)
library(irlba)

library(scchicFuncs)

# Explore some old analysis  ----------------------------------------------

inf.old <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K9me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"
load(inf.old, v=T)

out.lda <- out.lda[[1]]

topics.mat <- posterior(out.lda)$topics
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
         geom_point() + 
         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

lsi.out <- RunLSI(as.matrix(count.mat))

umap.out.lsi <- umap(lsi.out$u, config = jsettings)
dat.umap.long.lsi <- data.frame(cell = rownames(umap.out.lsi$layout), umap1 = umap.out.lsi$layout[, 1], umap2 = umap.out.lsi$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell))

ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2, color = experi)) + 
         geom_point() + 
         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load new mat ------------------------------------------------------------

inf.new <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K9me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
load(inf.new, v=T)

lsi.out.new <- RunLSI(as.matrix(count.mat))

umap.out.lsi.new <- umap(lsi.out.new$u, config = jsettings)
dat.umap.long.lsi.new <- data.frame(cell = rownames(umap.out.lsi.new$layout), umap1 = umap.out.lsi.new$layout[, 1], umap2 = umap.out.lsi.new$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell))

ggplot(dat.umap.long.lsi.new, aes(x = umap1, y = umap2, color = experi)) + 
         geom_point() + 
         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# add otal counts
count.sums <- data.frame(cell = colnames(count.mat), cellsize = colSums(count.mat))

dat.umap.long.lsi.new <- left_join(dat.umap.long.lsi.new, count.sums)


ggplot(dat.umap.long.lsi.new, aes(x = umap1, y = umap2, color = log10(cellsize))) + 
         geom_point() + 
         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

ggplot(dat.umap.long.lsi.new, aes(x = log10(cellsize)), fill = experi) + 
         geom_density() + 
         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()




# Load new LDA  -----------------------------------------------------------

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_H3K27me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K27me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"
load(inf, v=T)

tm.result <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

