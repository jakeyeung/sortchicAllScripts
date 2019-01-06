# Jake Yeung
# Date of Creation: 2018-12-28
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/lda_on_H3K27me3_bins.R
# Alexander found interesting Hox gene variability across cells using binning approach. Does this data show 
# information using alternative analyses?

rm(list=ls())

setwd("~/projects/scChiC")

library(dplyr)
library(ggplot2)
library(destiny)
library(umap)
library(princurve)
library(hash)

source("scripts/Rfunctions/FitFunctions.R")

# Load data ---------------------------------------------------------------

load("/private/tmp/lda_output/BM-H3K27me3.AvO_filt.Robj", v=T)
# load("/private/tmp/lda_output/lda_out.meanfilt.K-12.Robj", v=T)

# load("/private/tmp/lda_output/lda_outputs.meanfilt_1.merge_1000.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-12.Robj")
# load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)

load("/private/tmp/lda_output/lda_outputs.meanfilt_1.merge_1000_NoXYM.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-12.Robj", v=T)
load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)

count.mat <- count.dat$counts

# Plots -------------------------------------------------------------------


rownames(out.lda@gamma) <- out.lda@documents

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
jpeaks <- names(head(loads1, n = topn))[1:topn]
# jpeaks <- names(head(loads2, n = topn))[1:topn]
# jpeaks <- names(head(loadsqr, n = topn))[1:topn]
# jpeaks <- names(tail(loads1, n = topn))[1:topn]    # contains Hox in chr15 and chr11 
# (jpeaks <- grep("chr7:1035", names(loadsqr), value = TRUE))  # Hbb locus
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
rownames(dat.pca) <- rownames(out.lda.pca$rotation)

ggplot(dat.pca, aes(x = PC1, y = PC2, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC2, y = PC3, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

ggplot(dat.pca, aes(x = PC3, y = PC4, size = peak.counts, col = peak.counts)) + geom_point(alpha = 0.9)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))

# do UMAP

nn <- 20
jmetric <- 'pearson2'
jmindist <- 0.001
custom.settings <- umap.defaults
custom.settings$n_neighbors <- nn
custom.settings$metric <- jmetric
custom.settings$min_dist <- jmindist

dat.umap <- umap(out.lda@gamma, config = custom.settings)
# dat.umap <- umap(as.matrix(count.mat), config = custom.settings)
rownames(dat.umap$layout) <- rownames(out.lda@gamma)

rbPal <- colorRampPalette(c('black','yellow'))
jcol.rgb <- rbPal(20)[as.numeric(cut(jcol,breaks = 20))]

jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, col = jcol.rgb, main = jmain)


# do diffusion 
# 
# out.lda.dm <- DiffusionMap(out.lda@gamma, k = 12, n_eigs = 3)
# plot(out.lda.dm, 1:2, col = jcol)  # colors from beta matrix
# plot(out.lda.dm, 1:3, col = jcol)  # colors from beta matrix
# plot(out.lda.dm, 2:3, col = jcol)  # colors from beta matrix


# Run pseudotime and plot Hox genes ---------------------------------------

# jsub <- as.data.frame(out.lda.dm@eigenvectors) %>%
  # arrange(DC1)

jsub <- as.data.frame(dat.pca) %>%
  arrange(PC1)

rownames(jsub) <- out.lda@documents

# jsub <- dat.pca %>%
  # arrange(PC1)

plot(jsub[, 1], jsub[, 2])

# curvestart <- lm(DC2 ~ 1 + DC1 + I(DC1 ^ 2), jsub)
curvestart <- lm(PC2 ~ 1 + PC1 + I(PC1 ^ 2), jsub)

plot(jsub[, 1], jsub[, 2])
lines(jsub[, 1], curvestart$fitted.values)

startvals <- cbind(jsub[, 1], curvestart$fitted.values)

# plot(out.lda.dm@eigenvectors[, 2], curvestart$fitted.values)
  
proj <- principal_curve(as.matrix(jsub[, 1:2]), start = startvals)

dat.pca.proj <- as.data.frame(proj$s)
dat.pca.proj$lambda <- proj$lambda
dat.pca.proj <- dat.pca.proj %>%
  arrange(lambda)

plot(proj)
lines(proj)
points(proj)
whiskers(as.matrix(jsub), proj$s)

# # replot PCA 
# ggplot(jsub, aes(x = DC1, y = DC2)) + geom_point(alpha = 0.9)  +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings")) +
#   geom_line(data = dat.pca.proj, inherit.aes = FALSE, aes(x = DC1, y = DC2, alpha = lambda))

pseudohash <- hash(names(proj$lambda), proj$lambda)

# plot PCA with pseudotime as color

dat.pca$pseudo <- sapply(rownames(dat.pca), function(x) pseudohash[[x]])


ggplot(dat.pca, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.2)  +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings")) + 
  geom_line(data = dat.pca.proj, inherit.aes = FALSE, aes(x = PC1, y = PC2, alpha = lambda))


# Do counts regression ----------------------------------------------------

jpseudo <- sapply(colnames(count.mat), function(x) pseudohash[[x]])
jpseudo <- scale(jpseudo, center = TRUE, scale = TRUE)  # scale and center
jsize <- Matrix::colSums(count.mat)

# fit all data
print("Running fits...")

system.time(
  glm.fits <- apply(count.mat, MARGIN = 1, FUN = FitGlmRow, pseudo = jpseudo, size = jsize)
)
save(glm.fits, jpseudo, jsize, file = "/Users/yeung/projects/scChiC/outputs_R/fit_output/glm_fits_peak_calling_windows_withOffset.Robj")


# (jpeak.test <- grep("chr15:1029", rownames(count.mat), value = TRUE))
# jpeak.test.i <- which(rownames(count.mat) == jpeak.test[[9]])
# fit.test <- FitGlmRow(count.mat[jpeak.test.i, ], pseudo = jpseudo, size = jsize, returnobj = TRUE)
# 
# # plot fit
# # pred.x <- runif(100, min(jpseudo), max(jpseudo))
# # pred.y <- exp(predict(fit.test, newdata = data.frame(pseudo = pred.x, size = mean(jsize))))
# 
# pred.x <- jpseudo
# pred.y <- fits.test$fitted.values
# 
# jdat.real <- data.frame(x = jpseudo, y = count.mat[jpeak.test.i, ], type = "real")
# jdat.pred <- data.frame(x = pred.x, y = pred.y, type = "pred")
# 
# ggplot() + geom_point(data = jdat.real, aes(x = x, y = y), alpha=0.15) + 
#   geom_line(data = jdat.pred, aes(x = x, y = y))



# 
# ggplot() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings")) +
#   geom_line(data = dat.pca.proj, inherit.aes = FALSE, aes(x = PC1, y = PC2, alpha = lambda))
