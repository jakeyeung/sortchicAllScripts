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
library(MASS)
library(zinbwave)
library(DESeq2)
library(apeglm)
library(SummarizedExperiment)
library(GenomicRanges)

# Where do these peaks show up in the LDA model?  -------------------------

# load("/private/tmp/lda_output/BM-H3K27me3.AvO_filt.Robj", v=T)
load("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_tmp/BM-H3K27me3.AvO_filt.Robj", v=T)
load("/hpc/hub_oudenaarden/jyeung/data/scChiC/from_tmp/glm_fits_100kb_windows_withOffset.Robj", v=T)
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/zinb_deseq_output"

count.mat <- count.dat$counts

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

print("Running zinb")

zinb <- zinbwave(obj, K=0, BPPARAM=SerialParam(), epsilon=1e12)

save(zinb, file = file.path(outdir, "zinb_out.Robj"))

dds <- DESeqDataSet(zinb, design= ~pseudo)
# arguments as recommended from Van den Berge and Perraudeau
dds <- DESeq(dds, test="LRT", reduced=~1,
             sfType="poscounts", minmu=1e-6, minRep=Inf)

print("Running deseq2")

save(dds, file = file.path(outdir, "dds_out.Robj"))

zinb2 <- zinbwave(obj, K=2, BPPARAM=SerialParam(), epsilon=1e12)

save(zinb2, file = file.path(outdir, "/tmp/zinb_K2.Robj"))

print(Sys.time() - jstart)
