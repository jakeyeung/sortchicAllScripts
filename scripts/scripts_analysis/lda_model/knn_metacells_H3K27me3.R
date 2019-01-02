# Jake Yeung
# Date of Creation: 2019-01-01
# File: ~/Dropbox/scCHiC_figs/FIG4_BM/knn_metacells_H3K27me3.R
# Get colors for metacell, then link with LDA analysis


rm(list=ls())

jstart <- Sys.time() 

setwd("~/Dropbox/scChIC_figs/FIG4_BM")

library(hash)
library(ggplot2)
library(dplyr)
library(destiny)
library(Rtsne)
library(MASS)

# Functions ---------------------------------------------------------------

source("scripts/Rfunctions/FitFunctions.R")


# See MetaCell outputs ----------------------------------------------------



load("H3K27me3.datadir_mc_f.Rda")

Q<-object@mc_fp
mc_index<-object@mc
# mc_colors <- sapply(object@mc, function(x) ifelse(x %in% c(5, 7), "#A9A9A9", "#000000"))
mc_colors<-object@colors
colhash <- hash(seq(length(mc_colors)), mc_colors)
cellhash <- hash(names(mc_index), mc_index)
# mc_colors <- rep("#000000", 11)
# mc_colors[c(5, 7)] <- "#A9A9A9"

Q<-data.frame(Q)


load("H3K27me3_2d.datadir_2dproj.Rda")
x<-object@sc_x
y<-object@sc_y

# plot(x,y,pch=16,col=mc_colors[mc_index],cex=.75,axes=FALSE)

# rbPal <- colorRampPalette(c('black','yellow'))
# counts.total.col <- rbPal(100)[as.numeric(cut(log(counts.total),breaks = 100))]
# plot(x,y,pch=16,col=counts.total.col,cex=.75,axes=FALSE)




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
# jcol <- sapply(mc_index, function(x) colhash[[as.character(x)]])
# jcol <- sapply(mc_index, function(x) ifelse(x %in% c(5, 7), "#A9A9A9", "#000000"))


dat.pca <- as.data.frame(out.lda.pca$rotation) %>%
  mutate(counts.total = counts.total)  # if you want to see if there's bias with library size


dat.pca$peak.counts <- jcol
rownames(dat.pca) <- colnames(count.mat)

dat.pca$cluster <- sapply(rownames(dat.pca), function(x) cellhash[[x]])
dat.pca$is.5.7 <- sapply(dat.pca$cluster, function(x) ifelse(x %in% c(5, 7), TRUE, FALSE))
dat.pca$is.7 <- sapply(dat.pca$cluster, function(x) ifelse(x %in% c(7), TRUE, FALSE))
dat.pca$is.5 <- sapply(dat.pca$cluster, function(x) ifelse(x %in% c(5), TRUE, FALSE))

# ggplot(dat.pca, aes(x = PC1, y = PC2, col = as.character(cluster))) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))
# 
# ggplot(dat.pca, aes(x = PC1, y = PC2, col = counts.total)) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))
# 
# ggplot(dat.pca, aes(x = PC1, y = PC2, col = is.5.7)) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))
# 
# ggplot(dat.pca, aes(x = PC2, y = PC3, col = is.5.7)) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))
# 
# ggplot(dat.pca, aes(x = PC2, y = PC3, col = is.5)) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))


# Run principal curve -----------------------------------------------------

library(princurve)

proj <- principal_curve(as.matrix(dat.pca[, c("PC1", "PC2")]))

dat.pca.proj <- as.data.frame(proj$s)
dat.pca.proj$lambda <- proj$lambda
  

# ggplot(dat.pca, aes(x = PC1, y = PC2, col = is.5.7)) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings")) + 
#   geom_line(data = dat.pca.proj, inherit.aes = FALSE, aes(x = PC1, y = PC2, alpha = lambda))


# Plot peaks that show differences along pseudotime ------------------------

head(proj$lambda)

pseudohash <- hash(names(proj$lambda), proj$lambda)

# plot an example
# hoxc.peaks.i <- grep("chr15:103", rownames(count.mat), value = FALSE)
(hoxc.peaks <- grep("chr15:1029", rownames(count.mat), value = TRUE))
(hoxc.peaks.i <- grep("chr15:1029", rownames(count.mat), value = FALSE))
# hoxc.peaks.i <- sample(x = seq(nrow(count.mat)), size = 5, replace = FALSE)


# plot(x = sapply(colnames(count.mat), function(x) pseudohash[[x]]), colSums(count.mat[hoxc.peaks.i, ]), pch=20)

# plot(x = sapply(colnames(count.mat), function(x) pseudohash[[x]]), counts.total, pch=20, col = clstvec)

# fit using GLM


# dat <- data.frame(pseudo = sapply(colnames(count.mat), function(x) pseudohash[[x]]),
                  # counts = colSums(count.mat[hoxc.peaks.i, ]))
# m1 <- glm.nb(counts ~ 1 + pseudo, data = dat)
# m1.pois <- glm(counts ~ 1 + pseudo, data = dat, family = "poisson")
# m1.pois.null <- glm(counts ~ 1, data = dat, family = "poisson")
# 
# (est <- cbind(Estimate = coef(m1), confint(m1)))
# 
# m1.compare <- anova(m1.pois.null, m1.pois)
# # plot fits
# dat$phat <- m1.pois$fitted.values
# ggplot(dat) + geom_point(aes(x = pseudo, y = counts)) + geom_line(aes(x = pseudo, y = phat)) 

# Run fits for all genes

print("Running fits...")

jpseudo <- sapply(colnames(count.mat), function(x) pseudohash[[x]])
jpseudo <- scale(jpseudo, center = TRUE, scale = TRUE)  # scale and center
jsize <- Matrix::colSums(count.mat)
glm.fits <- apply(count.mat, MARGIN = 1, FUN = FitGlmRow, pseudo = jpseudo, size = jsize)
save(glm.fits, jpseudo, jsize, file = "/Users/yeung/projects/scChiC/outputs_R/fit_output/glm_fits_100kb_windows_withOffset.Robj")

print(Sys.time() - jstart)

# FitGlmRow(count.mat[1, ], pseudo = sapply(colnames(count.mat), function(x) pseudohash[[x]]))
# apply(count.mat[hoxc.peaks.i, ], MARGIN = 1, FitGlmRow, pseudo = jpseudo)





# Can we do pseudotime on this? -------------------------------------------

# library(umap)
# custom.settings <- umap.defaults
# custom.settings$n_neighbors <- 12
# custom.settings$metric <- "pearson"
# dat.umap <- umap(out.lda@gamma, config = custom.settings)
# rownames(dat.umap$layout) <- rownames(out.lda@gamma)
# clstvec <- sapply(rownames(out.lda@gamma), function(x) cellhash[[x]])
# is57 <- sapply(clstvec, function(x) ifelse(x %in% c(5, 7), "#A9A9A9", "#000000"))
# is7 <- sapply(clstvec, function(x) ifelse(x %in% c(7), "#A9A9A9", "#000000"))
# is5 <- sapply(clstvec, function(x) ifelse(x %in% c(5), "#A9A9A9", "#000000"))
# plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, col = is57)
# plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, col = clstvec)
# plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, col = is7)
# plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, col = is5)


# do UMAP on betas matrix
# betas.umap <- umap(t(exp(out.lda@beta)), config = custom.settings)
# plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20)

# count.pca <- prcomp(t(count.mat), center=TRUE, .scale=FALSE)

# Why is metacell 7 an outlier in UMAP analysis? --------------------------

# dat.tsne <- Rtsne(out.lda@gamma, dims = 2)
# plot(dat.tsne$Y, main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", pch = 20, col = clstvec)


# 
# 
# ggplot(dat.pca, aes(x = PC3, y = PC4, col = is.5.7)) + geom_point(alpha = 0.9)  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste0("Color by sum of top ", topn, " peaks from PC1 beta loadings"))
# 
# rownames(out.lda@gamma) <- out.lda@documents
# clstvec <- sapply(rownames(out.lda@gamma), function(x) cellhash[[x]])
# is57 <- sapply(clstvec, function(x) ifelse(x %in% c(5, 7), "#A9A9A9", "#000000"))
# out.lda.dm <- DiffusionMap(out.lda@gamma, k = 15)
# plot(out.lda.dm, 1:3, col = is57)  # colors from beta matrix
# plot(out.lda.dm, 1:2, col = clstvec)  # colors from beta matrix
# plot(out.lda.dm, 1:2, col = is7)  # colors from beta matrix
# plot(out.lda.dm, 2:3, col = is57) # colors from beta matrix
# 
# # color cells by Metacell 5 and 7 versus others
# 
# 
