# Jake Yeung
# Date of Creation: 2019-12-02
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_dedupfixed.R
# LDA dedup fixed

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(hash)
library(igraph)
library(umap)
library(Seurat)
library(scchicFuncs)

library(preprocessCore)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DESeq2)

library(ggrepel)

# Load data  --------------------------------------------------------------

jmark <- "H3K4me3"
# jsuff <- "UnenrichedXLinneg"
# jsuff <- "Unenriched"
# jsuff <- "AllMerged"
jsuff <- "Linneg"
jbin <- "TRUE"

# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_", 
#               jsuff, "_", jmark, "_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.binarize.", jbin, 
#               "/ldaOut.B6BM_", jsuff, "_", jmark, "_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.Robj")

inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_",  jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.binarize.", jbin, "/ldaOut.B6BM_", jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)  %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell))

# itereate umap
nns <- seq(from = 5, to = 100, by = 5)
mindists <- seq(from = 1, to = 5, by = 1)

pdf(paste0("/Users/yeung/data/scchic/pdfs/stemcell_linneg_analysis.red/AllMerged_test_umaps.pdf"), useDingbats = FALSE)



dev.off()
