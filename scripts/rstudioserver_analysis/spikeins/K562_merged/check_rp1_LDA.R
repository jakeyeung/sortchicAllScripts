# Jake Yeung
# Date of Creation: 2020-10-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_merged/check_rp1_LDA.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(Seurat)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "k36me3"

inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_RP1_ChICTAPS/lda_outputs.count_mat_from_seurat_stringent.k36me3.K-30.binarize.FALSE/ldaOut.count_mat_from_seurat_stringent.k36me3.K-30.Robj")
# inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_RP1_ChICTAPS/lda_outputs.count_mat_from_seurat.k36me3.K-30.binarize.FALSE/ldaOut.count_mat_from_seurat.k36me3.K-30.Robj")

load(inf.lda, v=T)

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result, jsep = "")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(data = dat.umap, mapping = aes(x = umap1, y = umap2)) + 
  geom_point()  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load indx ---------------------------------------------------------------


inf.indx <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/facs_index_data/RP1_FUCCI_2020-09-25/FACS_pseudotime_curve_plotted.with_hoescht.txt"
dat.indx <- as.data.frame(fread(inf.indx))
rownames(dat.indx) <- dat.indx$samp

dat.merge <- left_join(dat.umap, dat.indx, by = c("cell" = "samp"))


m0 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log2(X.405..460.50.Area))) + 
  geom_point()  + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("From LDA", jmark)


# Check raw input ---------------------------------------------------------

count.mat.norm <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")

pca.out <- prcomp(log2(t(count.mat.norm) + 1), retx = TRUE, center = TRUE)

dat.pca <- data.frame(cell = rownames(pca.out$x), pc1 = pca.out$x[, 1], pc2 = pca.out$x[, 2], stringsAsFactors = FALSE) %>%
  left_join(., dat.indx, by = c("cell" = "samp"))

dat.umap.frompca <- DoUmapAndLouvain(pca.out$x[, 1:50], jsettings) %>%
  left_join(., dat.indx, by = c("cell" = "samp"))

ggplot(dat.pca, aes(x = pc1, y = pc2, color = log2(X.405..460.50.Area))) + 
  geom_point() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() 

ggplot(dat.umap.frompca, aes(x = umap1, y = umap2, color = log2(X.405..460.50.Area))) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw()  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  


# Seurat ------------------------------------------------------------------

s.k562 <- CreateSeuratObject(counts = count.mat, project = paste0("k562_", jmark), min.cells = 10, min.features = 1000)
# s.k562 <- CreateSeuratObject(counts = dat.sub, project = paste0("k562_", jmark), min.cells = 0, min.features = 0)
s.k562 <- NormalizeData(s.k562, normalization.method = "LogNormalize", scale.factor = 10000)

s.k562 <- FindVariableFeatures(s.k562, selection.method = "vst", nfeatures = 25000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(s.k562), 10)

# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(s.k562)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# JFuncs::multiplot(plot1, plot2, cols = 1)

all.genes <- rownames(s.k562)
s.k562 <- ScaleData(s.k562, features = all.genes)
# s.k562 <- ScaleData(s.k562)

s.k562 <- RunPCA(s.k562, features = VariableFeatures(object = s.k562), npcs = 50)

DimPlot(s.k562, reduction = "pca")

s.k562 <- RunUMAP(s.k562, dims = 1:50)
DimPlot(s.k562, reduction = "umap")


umap.out <- data.frame(cell = rownames(s.k562@reductions$umap@cell.embeddings), s.k562@reductions$umap@cell.embeddings, stringsAsFactors = FALSE) %>%
  left_join(., dat.indx, by = c("cell" = "samp")) 

m1 <- ggplot(umap.out, aes(x = UMAP_1, y = UMAP_2, color = log2(X.405..460.50.Area))) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  ggtitle("From Seurat", jmark)

JFuncs::multiplot(m0, m1, cols = 2)


