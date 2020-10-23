# Jake Yeung
# Date of Creation: 2020-10-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_merged/compare_seurat_vs_LDA.RA
# 
# Check Seurat can give clusters? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(Seurat)
library(scchicFuncs)
library(JFuncs)

# Load tables -------------------------------------------------------------

# jmark <- "H3K4me1"
jmark <- "H3K27me3"
inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters/K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.countTable.binsize_100000.csv")

dat <- ReadMatSlideWinFormat(inf)

s.k562 <- CreateSeuratObject(counts = dat, project = paste0("k562_", jmark), min.cells = 3, min.features = 200)

s.k562 <- NormalizeData(s.k562, normalization.method = "LogNormalize", scale.factor = 10000)

s.k562 <- FindVariableFeatures(s.k562, selection.method = "vst", nfeatures = 200000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(s.k562), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(s.k562)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

JFuncs::multiplot(plot1, plot2, cols = 1)

all.genes <- rownames(s.k562)
s.k562 <- ScaleData(s.k562, features = all.genes)

s.k562 <- RunPCA(s.k562, features = VariableFeatures(object = s.k562))

DimPlot(s.k562, reduction = "pca")

s.k562 <- RunUMAP(s.k562, dims = 1:50)

DimPlot(s.k562, reduction = "umap")

dat.meta <- data.frame(cell = colnames(dat), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(indx = as.numeric(strsplit(cell, "_")[[1]][[2]]),
         row.indx = GetWellPosition(indx, platecols = 24, is.zero.base = FALSE)[[1]], 
         col.indx = GetWellPosition(indx, platecols = 24, is.zero.base = FALSE)[[2]],
         cellcycle = AddCellCycleLabel(colcoord = col.indx)) 

dat.umap <- data.frame(cell = colnames(s.k562), s.k562@reductions$umap@cell.embeddings, stringsAsFactors = FALSE) %>%
  left_join(., dat.meta)

ggplot(dat.umap, aes(x = UMAP_1, y = UMAP_2, color = cellcycle)) + 
  geom_point() + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Get index info ----------------------------------------------------------



# s.k562.sub <- subset(s.k562, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# 
# # Add FACS info  ----------------------------------------------------------
# 

inf.hoescht <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.hoesch_staining/hoesch_on_K562_plates.rds"
dat.hoescht <- readRDS(inf.hoescht)


dat.hoescht.all <- readRDS(inf.hoescht)