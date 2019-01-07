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
print(paste("Keeping", length(genes.keep), "genes"))

sumexp <- sumexp[genes.keep, ]

assayNames(sumexp)[1] <- "counts"

zinb <- zinbwave(sumexp, K = 2, epsilon=length(genes.keep))

W <- reducedDim(zinb)

