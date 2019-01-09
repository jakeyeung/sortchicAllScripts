# Jake Yeung
# Date of Creation: 2019-01-07
# File: ~/projects/scChiC/scripts/scripts_analysis/zinb_deseq/analyze_RNAseq_data.R
# Explore RNA-seq data

library(zinbwave)
library(data.table)

library(DESeq2)
library(apeglm)
library(SummarizedExperiment)
library(GenomicRanges)
library(ggplot2)
library(Seurat)
library(umap)
# library(velocyto.R)

# Load from Dropbox -------------------------------------------------------

inf <- "/Users/yeung/Dropbox/scCHiC_figs/BM-celseq2/PZ-BM-celseq2-countmatrix.dat"

dat <- data.table::fread(inf)


countmat <- as.matrix(subset(dat, select = -V1))
rownames(countmat) <- dat$V1

cellcounts <- colSums(countmat)
plot(density(log(cellcounts)))
abline(v = 5)
cells.remove <- cellcounts < exp(5)

print(dim(countmat))
countmat <- countmat[, !cells.remove]
print(dim(countmat))

coldat <- DataFrame(Batch = 1, Cell = colnames(countmat))

sumexp <- SummarizedExperiment(assays = list(counts = countmat),
                               colData = coldat)

xmean <- apply(sumexp@assays[[1]], 1, function(x) log10(mean(x)))
xvar <- apply(sumexp@assays[[1]], 1, function(x) log10(sqrt(var(x)) / mean(x)))

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(xmean, xvar, pch = 20, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1))
abline(a = 0, b = -0.5)


# Find variable genes -----------------------------------------------------

sobj <- CreateSeuratObject(countmat, project = "PZ-BM-Celseq2", normalization.method="LogNormalize")
sobj <- FindVariableGenes(sobj)


genes.keep <- sobj@var.genes
genes.erdr1 <- grep("Erdr1|Mid1", rownames(countmat), value = TRUE)
print(paste("Keeping", length(genes.keep), "genes"))

# plot original genes and show variable genes
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(xmean, xvar, pch = 20, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1))
text(xmean, xvar, labels = sapply(names(xmean), function(x) ifelse(x %in% genes.keep, x, "")))
# text(xmean, xvar, labels = sapply(names(xmean), function(x) ifelse(x %in% genes.erdr1, x, "")))
abline(a = 0, b = -0.5)


sumexp <- sumexp[genes.keep, ]

assayNames(sumexp)[1] <- "counts"

zinb <- zinbwave(sumexp, K = 2, epsilon=length(genes.keep))

W <- reducedDim(zinb)

nn <- 15
# jmetric <- 'pearson2'
# jmetric <- 'euclidean'
jmetric <- 'cosine'
jmindist <- 0.1
custom.settings <- umap.defaults
custom.settings$n_neighbors <- nn
custom.settings$metric <- jmetric
custom.settings$min_dist <- jmindist

W.umap <- umap(W, config = custom.settings)

dat.W <- data.frame(W, umap1 = W.umap$layout[, 1], umap2 = W.umap$layout[, 2])

ggplot(dat.W, aes(x = W1, y = W2)) + geom_point() + theme_classic()
ggplot(dat.W, aes(x = umap1, y = umap2)) + geom_point() + theme_classic()

# do UMAP directly on the matrix, use cosine

umap.out <- umap(assays(sumexp)[[1]], config = custom.settings)

dat.umap <- data.frame(umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

ggplot(dat.umap, aes(x = umap1, y = umap2)) + geom_point() + theme_classic()



# Find clusters -----------------------------------------------------------


