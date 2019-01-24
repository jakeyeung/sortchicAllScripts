# Jake Yeung
# Date of Creation: 2019-01-22
# File: ~/projects/scChiC/scripts/scripts_analysis/primetime/BM_H3K4me1_analysis.R
# Complete analysis of H3K4me1


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



# Load data ---------------------------------------------------------------


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

# 0.98 threshold 
# inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
# 0.96 threshold 
inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.0.96.Robj"))

inf.mc <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jchip, ".datadir_mc_f.Rda"))

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.GREAT))
assertthat::assert_that(file.exists(inf.mc))


load(inf, v=T)
load(inf.GREAT, v=T)

# normalize count.mat by total sum
count.mat <- sweep(count.mat, MARGIN = 2, STATS = Matrix::colSums(count.mat), FUN = "/") * 10^6

