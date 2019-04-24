# Jake Yeung
# Date of Creation: 2019-04-23
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_count_matrix.R
# Make count matrix

rm(list=ls())

library(Matrix)
library(dplyr)
library(ggplot2)
library(topicmodels)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)


source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")

# Correlate with peak expression ------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3")
jmarks.repress <- c("H3K27me3", "H3K9me3")
jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all

infs.nobin <- list("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
names(infs.nobin) <- c(jmarks, jmarks.repress)

out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs.nobin, SIMPLIFY = FALSE)
count.mat.lst <- lapply(out.objs.nobin, function(x) x$count.mat)

saveRDS(count.mat.lst, file = paste0("~/data/scchic/robjs/count_mat_lst_unnorm_cellmin_550_cellmax_500000.binarize.FALSE.", Sys.Date(), ".rds"))
