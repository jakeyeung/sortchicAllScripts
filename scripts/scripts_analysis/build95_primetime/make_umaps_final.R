# Jake Yeung
# Date of Creation: 2019-04-02
# File: ~/data/scchic/build95_primetime/make_umaps_final.R
# Make final UMAPs 

rm(list=ls())

library(GGally)
library(purrr)

library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

library(igraph)

library(GGally)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaData.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")


# Functions ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3")
jmarks.repress <- c("H3K27me3", "H3K9me3")
jmarks.all <- c(jmarks, jmarks.repress)
jmarks.sub <- jmarks

# jmarks.all <- jmarks.sub
names(jmarks.all) <- jmarks.all


barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))



# Load LDAs ---------------------------------------------------------------
# jmark <- "H3K4me1"
out.lda.new.lst <- list()
inf1 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-15_20_25_30_35.Robj")
load(inf1, v=T)
out.lda.new.lst[["H3K4me1"]] <- ChooseBestLDA(out.lda)
inf2 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-15_20_25_30_35.Robj")
load(inf2, v=T)
out.lda.new.lst[["H3K4me3"]] <- ChooseBestLDA(out.lda)
inf3 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj")
load(inf3, v=T)
out.lda.new.lst[["H3K27me3"]] <- ChooseBestLDA(out.lda)
inf4 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
load(inf4, v=T)
out.lda.new.lst[["H3K9me3"]] <- ChooseBestLDA(out.lda)

topics.mat.new.lst <- lapply(out.lda.new.lst, function(out.lda){
  return(posterior(out.lda)$topics)
})

infs.nobin <- list("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
names(infs.nobin) <- c(jmarks, jmarks.repress)

out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs.nobin, SIMPLIFY = FALSE)
count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

# waiting for repressive calculations to finish
out.lda.lst <- out.lda.new.lst[c("H3K4me1", "H3K4me3")]

infs <- list(inf1, inf2, inf3, inf4)
names(infs) <- c(jmarks, jmarks.repress)
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf, convert.chr20.21.to.X.Y = TRUE), jmarks.all, infs, SIMPLIFY = FALSE)

tm.result.lst <- lapply(out.lda.new.lst, function(x) posterior(x))
count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topic %*% tm.result$terms)))

annots.lst <- lapply(c(jmarks, jmarks.all), function(jmark) out.objs[[jmark]]$regions.annot)




# Get the UMAPs  ----------------------------------------------------------


