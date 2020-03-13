# Jake Yeung
# Date of Creation: 2020-03-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/make_merged_exprs_mat_for_MARA.R
# Make a new merged data (centered and uncentered) from new LDA

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)

# Load data  --------------------------------------------------------------

# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

dat.out <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains_marks_merged/lda_outputs.", jmark, "-BM_AllMerged.merged_by_clusters_no_NAs.K-30.binarize.FALSE/ldaOut.", jmark, "-BM_AllMerged.merged_by_clusters_no_NAs.K-30.Robj")
  load(inf, v=T)
  tm.result <- posterior(out.lda)
  dat.impute <- t(tm.result$topics %*% tm.result$terms)
  return(list(dat.impute = dat.impute, count.mat = count.mat, tm.result = tm.result))
})

dat.imputes <- lapply(dat.out, function(x) x$dat.impute)
dat.mats <- lapply(dat.out, function(x) x$dat.impute)

dat.imputes.merge <- do.call(cbind, dat.imputes)
dat.mats.merge <- do.call(cbind, dat.mats)


# Center and rorate -------------------------------------------------------

ntopics <- nrow(dat.out[[1]]$tm.result$topics)
# GLM loglink function for multinom is log( p / (1 - p) )
# V %*% t(U) on init matrix gives an estimate of p
# after estimate p / (1 - p), THEN remove mean and finally do SVD to get factors and loadings estimate
p <- dat.imputes.merge
logodds <- log(p / (1 - p))
# remove mean and SVD
logodds.centered <- t(scale(t(logodds), center = TRUE, scale = FALSE))
# # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
# logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
# U.init <- logodds.pca$x  # cells by k
# V.init <- logodds.pca$rotation  # genes by k, no need to transpose

# make nice
geneids <- sapply(rownames(logodds), function(x) strsplit(x, ";")[[1]][[1]])
logodds.out <- data.frame(Gene.ID = geneids, logodds, stringsAsFactors = FALSE)
logodds.centered.out <- data.frame(Gene.ID = geneids, logodds.centered, stringsAsFactors = FALSE)

# write logodds to output
outname <- paste0("ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30.txt")
outdir <- "/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/count_mats_peaks_norm_merged"
outdir2 <- "/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/count_mats_peaks_unnorm_LDA_merged"
outtxt <- file.path(outdir, outname)
outtxt2 <- file.path(outdir2, outname)

fwrite(logodds.centered.out, file = outtxt, sep = "\t", col.names = TRUE, row.names = FALSE)
fwrite(logodds.out, file = outtxt2, sep = "\t", col.names = TRUE, row.names = FALSE)
