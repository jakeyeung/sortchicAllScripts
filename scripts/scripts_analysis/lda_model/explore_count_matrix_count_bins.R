# Jake Yeung
# Date of Creation: 2019-02-21
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/explore_count_matrix_count_bins.R
# Count unique cuts

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

library(tidyr)
library(scales)

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

out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
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


# Show unique cuts per cell  ----------------------------------------------

jmark <- "H3K4me3"
jmark <- "H3K4me3"

dlong <- lapply(jmarks, function(jmark){
  cell.sizes <- Matrix::colSums(count.mat.lst[[jmark]])
  dlong <- data.frame(cell = names(cell.sizes), counts = cell.sizes, mark = jmark)
  return(dlong)
}) %>%
  bind_rows()

ggplot(dlong %>% filter(mark == "H3K27me3"), aes(x = log10(counts))) + geom_histogram(bins = 50) + facet_wrap(~mark, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) 

# dlong <- as.data.frame(as.matrix(count.mat.lst[[jmark]])) %>%
#   gather(cell, counts) %>%
#   group_by(cell) %>%
#   summarise(counts = sum(counts))
# 

