# Jake Yeung
# Date of Creation: 2019-02-16
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/explore_count_matrix.R
# Explore count matrix

rm(list=ls())

library(dplyr)
library(ggplot2)
library(JFuncs)
library(topicmodels)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(hash)
library(umap)
library(Matrix)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")

# Load LDA ----------------------------------------------------------------


jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
jcolvec <- c("blue", "gray80", "red")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
names(out.objs) <- jmarks
names(out.objs.nobin) <- jmarks

tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))

# get count matrices
count.mat.lst <- lapply(out.objs.nobin, function(x) x$count.mat)

# what's the distirbution of counts?
lapply(count.mat.lst, range)

# sparsity of matrix?

sparsity <- lapply(count.mat.lst, function(count.mat) 1 - nnzero(count.mat) / length(count.mat))
multicounts <- lapply(count.mat.lst, function(count.mat) length(which(count.mat > 1))  / length(count.mat))

plot(density(unlist(as.matrix(count.mat.lst[[3]]))))

# cells with large counts, what are they lonely ?

indx <- which(count.mat.lst[[1]] > 50, arr.ind = TRUE)

