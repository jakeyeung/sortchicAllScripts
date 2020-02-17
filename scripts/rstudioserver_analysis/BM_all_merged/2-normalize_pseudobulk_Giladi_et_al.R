# Jake Yeung
# Date of Creation: 2020-01-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/2-normalize_pseudobulk_Giladi_et_al.R
# Normalize pseudobulks

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(igraph)
library(umap)
library(DESeq2)
library(preprocessCore)

remove.celltype <- "C1qb"

inf <- "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.RData"
load(inf, v=T)
inf.annot <- "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/tier3_annotation.txt"
dat.annot <- as.data.frame(fread(inf.annot))

dat.sum <- dat.sum[, which(colnames(dat.sum) != remove.celltype)]

coldata <- data.frame(marker = colnames(dat.sum))
rownames(coldata) <- colnames(dat.sum)

dds <- DESeqDataSetFromMatrix(countData = dat.sum,
                              colData = coldata,
                              design= ~ 1)
dds <- DESeq(dds)

dat.sum.norm <- assay(vst(dds))

plot(density(dat.sum.norm))

plot(boxplot(dat.sum.norm))

# check Car1 indeed have highly expreessed Car1?
sort(dat.sum.norm[, "Car1"], decreasing = TRUE)[1:10]

dat.sum.norm["Car1", ]
dat.sum.norm["Prg2", ]
dat.sum.norm["Siglech", ]
dat.sum.norm["Ccl5", ]

# do normalization

dat.sum.norm.quantnorm <- preprocessCore::normalize.quantiles(dat.sum.norm, copy = TRUE)
rownames(dat.sum.norm.quantnorm) <- rownames(dat.sum.norm)
colnames(dat.sum.norm.quantnorm) <- colnames(dat.sum.norm)

plot(boxplot(dat.sum.norm.quantnorm))

save(dat.sum.norm, dat.sum.norm.quantnorm, ncells.vec, file = "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")

