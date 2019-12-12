# Jake Yeung
# Date of Creation: 2019-12-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/2-explore_Giladi_et_al.R
# Explore the scRNAseq dataset 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

Vectorize(plog2p <- function(p){
  return(ifelse(p == 0, 0, p * log2(p)))
}, vectorize.args = "p")

CalculateEntropy <- function(p, normalize.p = FALSE){
  if (normalize.p){
    p <- p / sum(p)
  }
  S <- -sum(plog2p(p))
  return(S)
}




# Functions ---------------------------------------------------------------


ReadGiladi <- function(inf, remove.first.col = TRUE){
  dat <- fread(inf) %>% 
    dplyr::rename(gene = V1)
  if (remove.first.col){
    dat$gene <- NULL
  }
  return(dat)
}


# Load marker genes  ------------------------------------------------------

inf.markers <- "/home/jyeung/data/from_cluster/public_data/Giladi_et_al_2018/41556_2018_121_MOESM4_ESM.markergenes.xlsx"
assertthat::assert_that(file.exists(inf.markers))

library(xlsx)

markers <- xlsx::read.xlsx(inf.markers, sheetIndex = 1)


# Load annotations from Giladi's email  -----------------------------------

inf.meta <- "~/hpc/scChiC/public_data/Giladi_et_al_2018/tier3_annotation.txt"
assertthat::assert_that(file.exists(inf.meta))

dat.meta <- fread(inf.meta)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.meta, aes(x = x, y = y, color = marker)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load data ---------------------------------------------------------------



indir <- "~/data/from_cluster/public_data/Giladi_et_al_2018"
inf.meta <- file.path(indir, "GSE92575_metadata.txt")
infs <- list.files(path = indir, pattern = "*.txt.gz", full.names = TRUE)
names(infs) <- sapply(infs, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

# get gene list
genes <- ReadGiladi(infs[[1]], remove.first.col = FALSE)$gene

dats <- lapply(infs, ReadGiladi, remove.first.col = TRUE) %>%
  bind_cols()

meta <- as.data.frame(fread(inf.meta, skip = 14))
rownames(meta) <- meta$well

cells <- meta$well

dats.filt <- dats[, ..cells]

dats.filt <- as.matrix(dats.filt)
rownames(dats.filt) <- genes

# split into tiers and then calculate entropy?

library(Seurat)


pbmc <- CreateSeuratObject(counts = dats.filt, project = "pbmc3k", min.cells = 3, min.features = 200, meta.data = meta)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", group.by = "tier") + scale_color_viridis_d()
DimPlot(pbmc, reduction = "umap", group.by = "tier", split.by = "tier", ncol = 5) + scale_color_viridis_d()


# Plot gene counts --------------------------------------------------------

stemcell.genes <- as.character(subset(markers, gene.module == "Stem genes")$gene)

jmodules <- c("Stem genes", "Monocyte", "Stage I neutrophil", "Stage II neutrophil")
marker.genes <- as.character(subset(markers, gene.module %in% jmodules)$gene)
marker.genes <- as.character(subset(markers, gene.module %in% jmodules)$gene)

# jgenes <- "Hlf"
# jgenes <- stemcell.genes[1:10]

jgenes <- "S100a8"
jgenes <- "Hlf"

# RidgePlot(pbmc, stemcell.genes, group.by = "tier")
# VlnPlot(pbmc, stemcell.genes, group.by = "tier")

FeaturePlot(pbmc, features = jgenes, order = TRUE)

# plot entropy in single cells


# x <- pbmc@assays$RNA@scale.data[, 1]

# CalculateEntropy(2^x, normalize.p = TRUE)

# S.vec <- apply(pbmc@assays$RNA@data, 2, function(jcell) CalculateEntropy(2^jcell, normalize.p = TRUE))
# S.vec <- apply(pbmc@assays$RNA@data, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
S.vec <- apply(pbmc@assays$RNA@counts, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
jmeta <- data.frame(S = S.vec)
rownames(jmeta) <- names(S.vec)

pbmc <- AddMetaData(object = pbmc, metadata = jmeta, col.name = "entropy")

FeaturePlot(pbmc, features = 'entropy') + scale_color_viridis_c(direction = 1)





