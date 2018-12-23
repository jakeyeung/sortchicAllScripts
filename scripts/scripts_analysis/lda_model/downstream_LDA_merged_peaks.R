# Jake Yeung
# Date of Creation: 2018-12-22
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_LDA_merged_peaks.R
# Analyze LDA output on merged peaks from 25kb distance
# If this works it serves as input to metacell so we can compare

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/Rfunctions/Aux.R")

# Load LDA output ---------------------------------------------------------

lda.path <- "outputs_R/lda_output/lda_outputs.meanfilt.merge_25000.Robj"
load(lda.path ,v=T)

outdir <- "/tmp/lda_output/lda_outputs.meanfilt.merge_25000"

# Plot clusters -----------------------------------------------------------


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
topn <- 5
jpeaks <- names(head(loads1, n = topn))[1:topn]  
(jpeaks <- grep("chr7:1035", names(loadsqr), value = TRUE))  # Hbb locus
jpeaks <- c()
jpeaks.i <- which(rownames(count.mat) %in% jpeaks)
jcol <- Matrix::colSums(count.mat[jpeaks.i, ])
counts.total <- Matrix::colSums(count.mat)
dat.pca <- as.data.frame(out.lda.pca$rotation) %>%
  mutate(counts.total = counts.total)  # if you want to see if there's bias with library size
dat.pca$peak.counts <- jcol
rownames(dat.pca) <- seq(ncol(count.mat))

ggplot(dat.pca, aes(x = PC1, y = PC2, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC2, y = PC3, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC3, y = PC4, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC4, y = PC5, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))


# Write count mat to file so we can do metacell ---------------------------

write.table(as.matrix(count.mat), 
            file = "/tmp/metacell_inputs/PZ-BM-H3K4me1.merged.25000kbmerge.mat", 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

