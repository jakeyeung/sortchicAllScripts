# Jake Yeung
# Date of Creation: 2019-12-02
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_dedupfixed.R
# LDA dedup fixed

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)


# FUNCTIONS ---------------------------------------------------------------


MergeCells <- function(mat, dat.split, jname, cname = "cell"){
  cells.keep <- as.character(dat.split[[jname]][[cname]])
  assertthat::assert_that(length(cells.keep) > 0)
  return(rowSums(mat[, cells.keep]))
}

# Settings ----------------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

# Load data  --------------------------------------------------------------

inf.linneg <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_Linneg_H3K4me3_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.binarize.TRUE/ldaOut.B6BM_Linneg_H3K4me3_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.Robj"
inf.wt <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.binarize.TRUE/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.Robj"

assertthat::assert_that(file.exists(inf.wt))
assertthat::assert_that(file.exists(inf.linneg))

load(inf.linneg, v=T)
count.mat.linneg <- count.mat.orig
lda.linneg <- out.lda
tm.result.linneg <- posterior(out.lda)

load(inf.wt, v=T)
count.mat.wt <- count.mat.orig
lda.wt <- out.lda
tm.result.wt <- posterior(out.lda)

# tm.result.lst <- list(wt = tm.result.wt$topics, linneg = tm.result.linneg$topics)
tm.result.lst <- list(wt = tm.result.wt, linneg = tm.result.linneg)

dat.umap.long.lst <- lapply(tm.result.lst, function(tm.result) DoUmapAndLouvain(tm.result$topics, jsettings = jsettings) %>% rowwise() %>% mutate(experi= ClipLast(as.character(cell))))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#0b1b7f", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01")
ggplot(dat.umap.long.lst[[1]], aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
ggplot(dat.umap.long.lst[[2]], aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

# check variance

dat.impute.log.lst <- lapply(tm.result.lst, function(x) log2(t(x$topics %*% x$terms)))

dat.var.lst <- lapply(dat.impute.log.lst, function(dat.impute.log) CalculateVarAll(dat.impute.log, jchromos) %>% mutate(experi = ClipLast(as.character(cell))))

# visualize all
dat.var.lst %>%
  bind_rows() %>%
  ggplot(., aes(x = cell.var.within.sum, fill = experi)) + geom_density(alpha = 0.5)


# visualize 
jnames <- names(tm.result.lst)
names(jnames) <- jnames
dat.umap.long.merge.lst <- lapply(jnames, function(x) left_join(dat.umap.long.lst[[x]], dat.var.lst[[x]]))


# pseudobulk louvains
ggplot(dat.umap.long.merge.lst[[1]], aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)

ggplot(dat.umap.long.merge.lst[[2]], aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)


# check the raw
# for WT: see WT versus 
louv.prog.wt <- "3"
louv.diff.wt <- "6"

louv.prog.linneg <- "3"
louv.diff.linneg <- "7"

jsplit.wt <- split(dat.umap.long.merge.lst$wt, f = dat.umap.long.merge.lst$wt$louvain)
jsplit.linneg <- split(dat.umap.long.merge.lst$linneg, f = dat.umap.long.merge.lst$linneg$louvain)

x.prog.wt <- MergeCells(count.mat.wt, jsplit.wt, louv.prog.wt)
x.diff.wt <- MergeCells(count.mat.wt, jsplit.wt, louv.diff.wt)

# check chromo 15
bins.filt <- grep("^chr15:", rownames(count.mat.wt), value = TRUE)
# sort by numeri
bins.filt.sorted <- bins.filt[order(sapply(bins.filt, GetStart), decreasing = FALSE)]

par(mfrow=c(2,1))
plot(log2(x.prog.wt[bins.filt.sorted] + 1), type = "l", main = "x.prog.wt")
plot(log2(x.diff.wt[bins.filt.sorted] + 1), type = "l", main = "x.diff.wt")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


# check Linneg
x.prog.linneg <- MergeCells(count.mat.linneg, jsplit.linneg, louv.prog.linneg)
x.diff.linneg <- MergeCells(count.mat.linneg, jsplit.linneg, louv.diff.linneg)

par(mfrow=c(2,1))
plot(log2(x.prog.linneg[bins.filt.sorted] + 1), type = "l", main = "x.prog.linneg")
plot(log2(x.diff.linneg[bins.filt.sorted] + 1), type = "l", main = "x.diff.linneg")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


# 
# topics.mat <- tm.result.wt$topics
# terms.mat <- tm.result.wt$terms
# umap.out <- umap(topics.mat, config = jsettings)
# dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)  %>%
#   rowwise() %>%
#   mutate(experi = ClipLast(cell))
#   
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # do variance
# 
# dat.impute.log <- log2(t(topics.mat %*% terms.mat))
# jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# 
# dat.var <- CalculateVarAll(dat.impute.log, jchromos)
# 
# dat.umap.long.merge <- left_join(dat.umap.long, dat.var)
# 
# ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_viridis_c(direction = -1) + facet_wrap(~experi)
# 
# ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # check whether we see anything in the raw data? 
# cells.linneg <- grepl(pattern = "Linneg", colnames(count.mat.orig))
# cells.wt <- grepl(pattern = "^B6-13W1", colnames(count.mat.orig))
# 
# count.mat.linneg <- count.mat.orig[, cells.linneg]
# count.mat.wt <- count.mat.orig[, cells.wt]




# # use entropy as a measurement??
# 
# 

# # Merge by louvains and plot  ---------------------------------------------
# 
# dat.umap.long.merge.louv <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.merge)
# 
# # visualize the variance??
# nnzero(count.mat.linneg) / length(count.mat.linneg)
# nnzero(count.mat.sc) / length(count.mat.sc)



