# Jake Yeung
# Date of Creation: 2019-01-22
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_analysis_H3K4me1_fewer_K.R
# User fewer K

rm(list=ls())

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(umap)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

jchip <- "H3K4me1"

jdist <- 1000L
jmean <- 1
jmin <- 100L
jmax <- 500000L
# binarize <- "TRUE";  jtops <- "5_7_10_12_15_20_25_30"
binarize <- "FALSE"; jtops <- "15_20_25_30_35"

jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
infbase <- basename(inf)
infbase <- strsplit(infbase, ".Robj")[[1]][[1]]

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# use 15 K
out.lda <- out.lda[[1]]

kchoose <- out.lda@k

tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics

# print(head(tm.result$topics))
# cluster the terms, but first get interesting terms
# top.regions <- unique(unlist(topic.regions))  # topic.regions defined by threshold, 98 or 96th percentile of top weights in each column of the betas matrix
# terms.mat <- t(tm.result$terms)[top.regions, ]


nn=5
nnterms=15
jmetric='euclidean' 
jmindist=0.1
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)


jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))


jtopic <- 14
cell.assign <- apply(topics.mat, 1, which.max)

weights.in <- as.integer(cell.assign == jtopic)
weights.innot <- as.integer(cell.assign != jtopic)

col.binary <- weights.in + 1
# see where jtopic is in the map
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), 
     col = col.binary, asp = 0.75)