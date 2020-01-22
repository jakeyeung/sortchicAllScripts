# Jake Yeung
# Date of Creation: 2019-12-02
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_variance_by_pseudobulk.R
# Pseodbulk variance check

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)



# Settings ----------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load data  --------------------------------------------------------------

# load StemCells
jmark <- "H3K4me3"
jbin <- "TRUE"
jsuff <- "StemCells"

inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_",  jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.binarize.", jbin, "/ldaOut.B6BM_", jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)
count.mat.sc <- count.mat.orig

umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
dat.umap.long.louv <- DoLouvain()

# load lineage neg
jsuff <- "Linneg"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_",  jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.binarize.", jbin, "/ldaOut.B6BM_", jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)
count.mat.linneg <- count.mat.orig

# load unenriched
inf.unenriched <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.binarize.TRUE/ldaOut.B6BM_Unenriched_H3K4me3_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.Robj"
assertthat::assert_that(file.exists(inf.unenriched))
load(inf.unenriched, v=T)
count.mat.unenriched <- count.mat.orig
 
 

# Check sparsity  ---------------------------------------------------------

lapply(list(count.mat.unenriched, count.mat.linneg, count.mat.sc), function(x) nnzero(x) / length(x))



