# Jake Yeung
# Date of Creation: 2020-10-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/1-run_LDA_downstream_on_each_mark_and_dbl.filter_cells_downstream.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

# jmarjmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmarks <- c("H3K4me1", "H3K9me3", "H3K4me1xH3K9me3"); names(jmarks) <- jmarks
# hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load  -------------------------------------------------------------------

infdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_LDA_downstream.filt/LDA_downstream_objects.2020-10-21.again.RData"
outdir.rds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.single_match_dbl/cellfilt_binfilt"
dir.create(outdir.rds)

load(infdata, v=T)

# Write the topics --------------------------------------------------------

jtopics.k9me3 <- c("topic28", "topic7", "topic15", "topic20", "topic4")
names(jtopics.k9me3) <- jtopics.k9me3

jtopics.k4me1 <- c("topic16", "topic29", "topic24", "topic7", "topic18", "topic1", "topic8", "topic10")
names(jtopics.k4me1) <- jtopics.k4me1

# H3K9me3 topics
topn <- 500

jtopics.k9me3.vec <- lapply(jtopics.k9me3, function(jtopic){
  terms.tmp <- sort(tm.result.lst$H3K9me3$terms[jtopic, ], decreasing = TRUE)
  terms.tmp.filt <- terms.tmp[1:topn]
  return(names(terms.tmp.filt))
}) %>%
  unlist()
  
jtopics.k4me1.vec <- lapply(jtopics.k4me1, function(jtopic){
  terms.tmp <- sort(tm.result.lst$H3K4me1$terms[jtopic, ], decreasing = TRUE)
  terms.tmp.filt <- terms.tmp[1:topn]
  return(names(terms.tmp.filt))
}) %>%
  unlist()

bins.keep <- unique(c(jtopics.k9me3.vec, jtopics.k4me1.vec))


# Filter low var cells  ---------------------------------------------------


dat.merge.lst <- split(dat.merge, f = dat.merge$mark)
varfilt.lst <- list("H3K4me1" = 3, "H3K9me3" = 2.2, "H3K4me1xH3K9me3" = 0.8)

dat.filt.lst <- lapply(jmarks, function(jmark){
  varfilt <- varfilt.lst[[jmark]]
  dat.sub <- subset(dat.merge.lst[[jmark]], cell.var.within.sum.norm >= varfilt)
})

cells.filt.lst <- lapply(dat.filt.lst, function(jdat){
  jdat$cell
})

for (jmark in jmarks){
  print(jmark)
  print(paste0(nrow(dat.filt.lst[[jmark]]), "/", nrow(dat.merge.lst[[jmark]])))
  print(signif(nrow(dat.filt.lst[[jmark]]) / nrow(dat.merge.lst[[jmark]]), digits = 2))
}


# Write new count table ---------------------------------------------------


count.mat.filt.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  cells.filt <- cells.filt.lst[[jmark]]
  print(dim(count.mat.lst[[jmark]]))
  jmat <- count.mat.lst[[jmark]][bins.keep, cells.filt]
  print(dim(jmat))
  return(jmat)
})


# Write new output --------------------------------------------------------


for (jmark in jmarks){
  print(jmark)
  outname <- paste0("count_mat.", jmark, ".match_dbl.cellfilt.binfilt.rds")
  outf.rds <- file.path(outdir.rds, outname)
  saveRDS(count.mat.filt.lst[[jmark]], file = outf.rds)
}

# 
# library(irlba)
# 
# lsi.out <- RunLSI(as.matrix(count.mat.filt.lst$H3K4me1))
# 
# dat.umap.lsi <- DoUmapAndLouvain(lsi.out$u, jsettings)
# 
# ggplot(dat.umap.lsi, aes(x = umap1, y = umap2, color = louvain)) + geom_point() 
# 

