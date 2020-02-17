# Jake Yeung
# Date of Creation: 2019-12-01
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_outputs_dedupbug_vs_fixed.R
# Check LDA output

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


# Load dat ----------------------------------------------------------------

# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_nodedup/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-29.K-30.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_nodedup/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-29.K-30.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.binarize.TRUE/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.binarize.TRUE/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.Robj"
load(inf, v=T)

# LSI on original count mat
if (class(count.mat.orig) == "Matrix"){
  lsi.out <- RunLSI(as.matrix(count.mat.orig))
} else {
  lsi.out <- RunLSI(as.matrix(count.mat))
}

topics.mat <- posterior(out.lda)$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out <- umap(topics.mat, config = jsettings)
umap.out.lsi <- umap(lsi.out$u, config = jsettings)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
dat.umap.long.lsi <- data.frame(cell = rownames(umap.out.lsi$layout), umap1 = umap.out.lsi$layout[, 1], umap2 = umap.out.lsi$layout[, 2], stringsAsFactors = FALSE)

m.umap.lda <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point()  + ggtitle("Umap LDA")
m.umap.lda.lsi <- ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2)) + geom_point()  + ggtitle("Umap LSI")

print(m.umap.lda)
print(m.umap.lda.lsi)


# filter for good cells
# load old LDA
inf.dat.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
assertthat::assert_that(file.exists(inf.dat.stringent))
load(inf.dat.stringent, v=T)

cells.keep <- colnames(out.objs$count.mat)
bins.keep <- paste("chr", rownames(out.objs$count.mat), sep = "")

mat.sub <- count.mat.orig[, cells.keep]  # cells.keep is critical
lsi.out.filt <- RunLSI(as.matrix(mat.sub))

umap.out.lsi.filt <- umap(lsi.out.filt$u, config = jsettings)
dat.umap.long.lsi.filt <- data.frame(cell = rownames(umap.out.lsi.filt$layout), umap1 = umap.out.lsi.filt$layout[, 1], umap2 = umap.out.lsi.filt$layout[, 2], stringsAsFactors = FALSE)
m.umap.lda.lsi.filt <- ggplot(dat.umap.long.lsi.filt, aes(x = umap1, y = umap2)) + geom_point()  + ggtitle("Umap LSI, bin and cell filt")

print(m.umap.lda.lsi.filt)

# check density of good cells versus all cells

plot(density(colSums(count.mat.orig)))
plot(density(colSums(count.mat.orig[, cells.keep])))


