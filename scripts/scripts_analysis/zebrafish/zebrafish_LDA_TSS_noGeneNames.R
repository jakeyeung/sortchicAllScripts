# Jake Yeung
# Date of Creation: 2019-11-09
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_TSS.R
# Zebrafish LDA TSS

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

GetGeneAnnotsHash <- function(inf.annot){
  dat.annot <- data.table::fread(inf.annot, col.names = c("chromo", "start", "end", "bname"))
  # add chr
  dat.annot$chromo <- paste("chr", dat.annot$chromo, sep = "")
  rnames.old <- paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":")
  rnames.new <- dat.annot$bname
  annots.hash <- hash::hash(rnames.old, rnames.new)
}

AddGeneNameToRows <- function(mat, annots.hash){
  # mat rownmaes got stripped of gene names, add them back
  rnames.old <- rownames(mat)
  rnames.new <- sapply(rnames.old, function(x) annots.hash[[x]])
  rownames(mat) <- rnames.new
  return(mat)
}

# Load data ---------------------------------------------------------------

jmark <- "H3K4me1"
winsize <- 50000L

inf.annot <- paste0("/Users/yeung/data/scchic/tables/gene_tss.winsize_", winsize, ".species_drerio.nochr.bed")

# init
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow/lda_outputs.PZ-ChIC-ZFWKM-", jmark, 
              ".winsize_", winsize, ".merged.K-30.binarize.FALSE/lda_out_meanfilt.PZ-ChIC-ZFWKM-", jmark, 
              ".winsize_", winsize, ".merged.K-30.Robj")
load(inf, v=T)

out.lda <- out.lda[[1]]

print(head(out.lda@terms))

annots.hash <- GetGeneAnnotsHash(inf.annot)

count.mat <- AddGeneNameToRows(count.mat, annots.hash)

# out.lda@terms <- sapply(out.lda@terms, function(x) annots.hash)

# add gene name to the coordinates (got lost in mat to sparse mat pipeline)
topics.mat <- posterior(out.lda)$topics
terms.mat <- posterior(out.lda)$terms

colnames(terms.mat) <- sapply(colnames(terms.mat), function(x) annots.hash[[x]])

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R






# Find celltypes ----------------------------------------------------------




