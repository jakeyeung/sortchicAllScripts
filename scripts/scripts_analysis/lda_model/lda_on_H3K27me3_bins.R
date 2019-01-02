# Jake Yeung
# Date of Creation: 2018-12-28
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/lda_on_H3K27me3_bins.R
# Alexander found interesting Hox gene variability across cells using binning approach. Does this data show 
# information using alternative analyses?

rm(list=ls())


library(dplyr)
library(ggplot2)
library(destiny)

# Load data ---------------------------------------------------------------

# load("/private/tmp/lda_output/BM-H3K27me3.AvO_filt.Robj", v=T)
# load("/private/tmp/lda_output/lda_out.meanfilt.K-12.Robj", v=T)

# load("/private/tmp/lda_output/lda_outputs.meanfilt_1.merge_1000.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-12.Robj")
# load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)

load("/private/tmp/lda_output/lda_outputs.meanfilt_1.merge_1000_NoXYM.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-12.Robj", v=T)
load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)

count.mat <- count.dat$counts

# Plots -------------------------------------------------------------------


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
# jpeaks <- names(head(loads1, n = topn))[1:topn]
# jpeaks <- names(head(loads2, n = topn))[1:topn]
# jpeaks <- names(head(loadsqr, n = topn))[1:topn]
# jpeaks <- names(tail(loads1, n = topn))[1:topn]    # contains Hox in chr15 and chr11 
(jpeaks <- grep("chr7:1035", names(loadsqr), value = TRUE))  # Hbb locus
# jpeaks <- c()

# Hox cluster
# jpeaks <- grep("chr6:5213|chr2:7476|chr15:1029|chr11:9628", names(loads1), value = TRUE)

# chr2:74,760,980-74,767,142, chr15:102,964,796-102,973,897

print(jpeaks)

jpeaks.i <- which(rownames(count.mat) %in% jpeaks)
jcol <- Matrix::colSums(count.mat[jpeaks.i, ])
counts.total <- Matrix::colSums(count.mat)

dat.pca <- as.data.frame(out.lda.pca$rotation) %>%
  mutate(counts.total = counts.total)  # if you want to see if there's bias with library size

dat.pca$peak.counts <- jcol
rownames(dat.pca) <- seq(ncol(count.mat))

ggplot(dat.pca, aes(x = PC1, y = PC2, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC2, y = PC3, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC3, y = PC4, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

# do diffusion 

out.lda.dm <- DiffusionMap(out.lda@gamma, k = 15)
plot(out.lda.dm, 1:3, col = jcol)  # colors from beta matrix
plot(out.lda.dm, 1:2, col = jcol)  # colors from beta matrix
plot(out.lda.dm, 2:3, col = jcol)  # colors from beta matrix

# Can we get Hox clusters here? -------------------------------------------


