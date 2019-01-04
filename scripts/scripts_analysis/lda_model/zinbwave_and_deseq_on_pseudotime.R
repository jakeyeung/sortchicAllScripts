# Jake Yeung
# Date of Creation: 2019-01-01
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/zinbwave_and_deseq_on_pseudotime.R
# Use LDA to identify pseudotime, then run GLM model to find interesting genes


rm(list=ls())

jstart <- Sys.time() 

setwd("~/projects/scChiC")

library(hash)
library(ggplot2)
library(dplyr)
library(destiny)
library(Rtsne)
library(MASS)
library(princurve)
library(zinbwave)
library(DESeq2)
library(apeglm)
library(SummarizedExperiment)
library(GenomicRanges)

# Where do these peaks show up in the LDA model?  -------------------------

load("/private/tmp/lda_output/BM-H3K27me3.AvO_filt.Robj", v=T)
load("/private/tmp/lda_output/lda_out.meanfilt.K-12.Robj", v=T)

rownames(out.lda@gamma) <- out.lda@documents

count.mat <- count.dat$counts

counts.total <- Matrix::colSums(count.mat)

cells.keep <- out.lda@documents
peaks.keep <- out.lda@terms

cells.keep.i <- which(colnames(count.mat) %in% cells.keep)
peaks.keep.i <- which(rownames(count.mat) %in% peaks.keep)

count.mat <- count.mat[peaks.keep.i, cells.keep.i]

out.lda.pca <- prcomp(t(out.lda@gamma), center = TRUE, scale. = FALSE)
# betas are topic-to-peak weights. They are in natural log scale so I exponentiate. 
# Sum across peaks for each topic = 1
out.lda.betas.pca <- prcomp(exp(out.lda@beta), center = TRUE, scale. = FALSE) 

rownames(out.lda.betas.pca$rotation) <- out.lda@terms
rownames(out.lda.betas.pca$x) <- seq(nrow(out.lda.betas.pca$x))

# Plot using the gamma matrix using PCA or your favorite dim-reduc method. 
# You can also do PCA on the beta matrix, which clusters peaks
loads1 <- sort(out.lda.betas.pca$rotation[, 1], decreasing = TRUE)
loads2 <- sort(out.lda.betas.pca$rotation[, 2], decreasing = TRUE)
loads3 <- sort(out.lda.betas.pca$rotation[, 3], decreasing = TRUE)
loadsqr <- sort(sqrt(out.lda.betas.pca$rotation[, 1] ^ 2 + out.lda.betas.pca$rotation[, 2] ^ 2), decreasing = TRUE)


# color by top  peaks from PC1 of betas
topn <- 50
(jpeaks <- grep("chr7:1035", names(loadsqr), value = TRUE))  # Hbb locus
# jpeaks <- c()


print(jpeaks)

jpeaks.i <- which(rownames(count.mat) %in% jpeaks)
jcol <- Matrix::colSums(count.mat[jpeaks.i, ])

dat.pca <- as.data.frame(out.lda.pca$rotation) %>%
  mutate(counts.total = counts.total)  # if you want to see if there's bias with library size

dat.pca$peak.counts <- jcol
rownames(dat.pca) <- colnames(count.mat)

proj <- principal_curve(as.matrix(dat.pca[, c("PC1", "PC2")]))

dat.pca.proj <- as.data.frame(proj$s)
dat.pca.proj$lambda <- proj$lambda
  

# Plot peaks that show differences along pseudotime ------------------------

head(proj$lambda)

pseudohash <- hash(names(proj$lambda), proj$lambda)

jpseudo <- sapply(colnames(count.mat), function(x) pseudohash[[x]])

# plot an example
# hoxc.peaks.i <- grep("chr15:103", rownames(count.mat), value = FALSE)
(hoxc.peaks <- grep("chr15:1029", rownames(count.mat), value = TRUE))
(hoxc.peaks.i <- grep("chr15:1029", rownames(count.mat), value = FALSE))
# hoxc.peaks.i <- sample(x = seq(nrow(count.mat)), size = 5, replace = FALSE)


# Fit using zinb and deseq2 -----------------------------------------------

jseqnames <- sapply(rownames(count.mat), function(x) strsplit(x, ":")[[1]][[1]], USE.NAMES = FALSE)
jstartsends <- sapply(rownames(count.mat), function(x) strsplit(x, ":")[[1]][[2]], USE.NAMES = FALSE)
jstarts <- as.numeric(sapply(jstartsends, function(x) strsplit(x, "-")[[1]][[1]]))
jends <- as.numeric(sapply(jstartsends, function(x) strsplit(x, "-")[[1]][[2]]))

rowRanges <- GRanges(seqnames = jseqnames, ranges = IRanges(jstarts, jends))

colData <- data.frame(rname = colnames(count.mat), pseudo = jpseudo, row.names = NULL)

jcount.mat <- count.mat
rownames(jcount.mat) <- NULL
obj <- SummarizedExperiment(assays = list(counts = as.matrix(jcount.mat)),
                            rowRanges = rowRanges, colData = colData)

jstart <- Sys.time()

zinb <- zinbwave(obj, K=0, BPPARAM=SerialParam(), epsilon=1e12)

dds <- DESeqDataSet(zinb, design= ~pseudo)
# arguments as recommended from Van den Berge and Perraudeau
dds <- DESeq(dds, test="LRT", reduced=~1,
             sfType="poscounts", minmu=1e-6, minRep=Inf)

save(zinb, dds, file = "/Users/yeung/projects/scChiC/outputs_R/fit_output/zinb_dds_out.Robj")

zinb2 <- zinbwave(obj, K=2, BPPARAM=SerialParam(), epsilon=1e12)

save(zinb2, file = "/Users/yeung/projects/scChiC/outputs_R/fit_output/zinb_K2.Robj")

print(Sys.time() - jstart)
