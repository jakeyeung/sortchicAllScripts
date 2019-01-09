# Jake Yeung
# run_zinbwave.R
# Load count matrix, filter variable genes, run zinbwave, save object 
# 2019-01-07

jstart <- Sys.time()

library(zinbwave)
library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(ggplot2)
library(Seurat)
library(JFuncs)
library(BiocParallel)

args <- commandArgs(trailingOnly=TRUE)

inf <- args[[1]]
outf <- args[[2]]
k <- JFuncs::StrToNumeric(args[[3]])
ncores <- JFuncs::StrToNumeric(args[[4]])

dat <- data.table::fread(inf)

countmat <- as.matrix(subset(dat, select = -V1))
rownames(countmat) <- dat$V1

# remove bad cells?
cellcounts <- colSums(countmat)
plot(density(log(cellcounts)))
abline(v = 5)
cells.remove <- cellcounts < exp(5)

print(paste("Removing", length(cells.remove), "cells"))
print(dim(countmat))
countmat <- countmat[, !cells.remove]
print(dim(countmat))

coldat <- DataFrame(Batch = 1, Cell = colnames(countmat))

sumexp <- SummarizedExperiment(assays = list(counts = countmat),
                               colData = coldat)

xmean <- apply(sumexp@assays[[1]], 1, function(x) log10(mean(x)))
xvar <- apply(sumexp@assays[[1]], 1, function(x) log10(sqrt(var(x)) / mean(x)))

# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# plot(xmean, xvar, pch = 20, col = rgb(red = 1, green = 0, blue = 0, alpha = 0.1))
# abline(a = 0, b = -0.5)


# Find variable genes -----------------------------------------------------

sobj <- CreateSeuratObject(countmat, project = "PZ-BM-Celseq2", normalization.method="LogNormalize")
sobj <- FindVariableGenes(sobj)

genes.keep <- sobj@var.genes
print(paste("Keeping", length(genes.keep), "genes"))

sumexp <- sumexp[genes.keep, ]

assayNames(sumexp)[1] <- "counts"

# zinb <- zinbwave(sumexp, K = 2, epsilon=length(genes.keep))
zinb <- zinbFit(sumexp, K=k, epsilon=length(genes.keep), BPPARAM=MulticoreParam(ncores))

save(zinb, sobj, file = outf)
print(Sys.time() - jstart)

