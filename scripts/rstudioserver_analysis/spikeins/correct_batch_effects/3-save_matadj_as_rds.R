# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/1-correct_batch_effects_from_peaks_from_LDA.R
# description
rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(JFuncs)

library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DescTools)

library(topicmodels)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_batch_correction_output"
outdir <- indir

# make matrix by mvoing name to rowname, then removing rname colname
mat.adj.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir, paste0("mat_adj.", jmark, ".from_LDA.binskeep_1000.RData"))
  load(inf.tmp, v=T)  # mat.adj
  mat.adj <- as.data.frame(mat.adj)
  rownames(mat.adj) <- mat.adj$rname
  mat.adj$rname <- NULL
  return(mat.adj)
})


# Write outputs -----------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(outdir, paste0("mat_adj.", jmark, ".from_LDA.binskeep_1000.rds"))
  assertthat::assert_that(!file.exists(outrds))
  print(dim(mat.adj.lst[[jmark]]))
  print(mat.adj.lst[[jmark]][1:5, 1:5])
  saveRDS(mat.adj.lst[[jmark]], outrds)
}
