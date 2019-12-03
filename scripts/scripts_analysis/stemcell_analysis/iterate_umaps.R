# Jake Yeung
# Date of Creation: 2019-12-03
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/iterate_umaps.R
# 

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

library(topicmodels)

# Load data  --------------------------------------------------------------

jmark <- "H3K4me3"
# jsuff <- "UnenrichedXLinneg"
# jsuff <- "Unenriched"
jsuff <- "AllMerged"
# jsuff <- "Linneg"
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
  mutate(experi = ClipLast(cell)) %>%
  mutate(experi = gsub("-G2", "", experi))

# itereate umap
nns <- seq(from = 5, to = 100, by = 5)
mindists <- 10^-seq(from = 1, to = 5, by = 1)

pdf(paste0("/Users/yeung/data/scchic/pdfs/stemcell_linneg_analysis.redo/AllMerged_test_umaps.pdf"), useDingbats = FALSE)


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#ff9f7d", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
for (nn in nns){
  for (mindist in mindists){
	print(mindist)
    jsettings$min_dist <- mindist
    jsettings$n_neighbors <- nn
    umap.out <- umap(topics.mat, config = jsettings)
    dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)  %>%
      rowwise() %>%
      mutate(experi = ClipLast(cell))
    m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_manual(values = cbPalette) + ggtitle(paste("Mindist:", mindist, "NN:", nn))
    print(m)
  }
}
dev.off()
