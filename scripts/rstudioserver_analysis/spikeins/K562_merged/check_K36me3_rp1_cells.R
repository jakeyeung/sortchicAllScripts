# Jake Yeung
# Date of Creation: 2020-10-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_merged/check_K36me3_rp1_cells.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)
library(scchicFuncs)

writeoutputs <- FALSE

hubpath <- "/home/jyeung/hub_oudenaarden"
# jmark <- "k27me3"
# jmark <- "k9me3"
jmark <- "k36me3"

jmarks <- c("k36me3", "k27me3", "k9me3"); names(jmarks) <- jmarks


# inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/merged_by_mark/countTablesAndRZr1only_ByChromo.NewFilters/CG-ChIC-TAPS-RPE-exp1-k36me3.tagged.countTable.binsize_100000.csv"))
inf <- file.path(hubpath, paste0("jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/merged_by_mark/countTablesAndRZr1only_ByChromo.NewFilters/CG-ChIC-TAPS-RPE-exp1-", jmark, ".tagged.countTable.binsize_100000.csv"))
# inf <- file.path(hubpath, paste0("avo/RPE1-ChIC/raw_demultiplexed/CG-qcChIC-RPE-exp1-", jmark, "-pl1/countmatrix"))
assertthat::assert_that(file.exists(inf))

dat <- ReadMatSlideWinFormat(inf)

# Filter bad cells  -------------------------------------------------------

plot(density(log10(colSums(dat))))

cell.sizes <- data.frame(cell = colnames(dat), ncuts = colSums(dat), stringsAsFactors = FALSE)
rownames(cell.sizes) <- colnames(dat)

counts.min <- 1000
counts.max <- 100000

cells.keep <- subset(cell.sizes, ncuts >= counts.min & ncuts <= counts.max)$cell

dat.sub <- dat[, cells.keep]

# Load in seurat  ---------------------------------------------------------


s.k562 <- CreateSeuratObject(counts = dat.sub, project = paste0("k562_", jmark), min.cells = 10, min.features = 1000)
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

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/facs_index_data/RP1_FUCCI_2020-09-25/FACS_pseudotime_curve_plotted.txt"
dat.indx <- as.data.frame(fread(inf))
rownames(dat.indx) <- dat.indx$samp


s.k562 <- AddMetaData(object = s.k562, metadata = dat.indx)
s.k562 <- AddMetaData(object = s.k562, metadata = cell.sizes)

umap.out <- data.frame(cell = rownames(s.k562@reductions$umap@cell.embeddings), s.k562@reductions$umap@cell.embeddings, stringsAsFactors = FALSE) %>%
  left_join(., dat.indx, by = c("cell" = "samp")) %>%
  left_join(., cell.sizes, by = c("cell"))

pca.out <- data.frame(cell = rownames(s.k562@reductions$pca@cell.embeddings), s.k562@reductions$pca@cell.embeddings, stringsAsFactors = FALSE) %>%
  left_join(., dat.indx, by = c("cell" = "samp")) %>%
  left_join(., cell.sizes, by = c("cell"))

ggplot(umap.out, aes(x = UMAP_1, y = UMAP_2, color = lambda)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(umap.out, aes(x = UMAP_1, y = UMAP_2, color = log10(ncuts))) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(pca.out, aes(x = PC_1, y = PC_2, color = lambda)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# Get input of PCA  -------------------------------------------------------

mat.input <- s.k562@assays$RNA@scale.data
cells.keep <- colnames(mat.input)

cells.keep.stringent <- subset(umap.out, UMAP_1 < 5)$cell
genes.keep <- rownames(mat.input)
pca.out <- prcomp(t(mat.input))

dat.pca <- data.frame(cell = rownames(pca.out$x), pc1 = pca.out$x[, 1], pc2 = pca.out$x[, 2], stringsAsFactors = FALSE)

ggplot(dat.pca, aes(x = pc1, y = pc2)) + 
  geom_point() + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# 
# # Write output for LDA ----------------------------------------------------
# 
# dat.sub.forLDA <- dat.sub[genes.keep, cells.keep]
# dat.sub.forLDA.stringent <- dat.sub[genes.keep, cells.keep.stringent]
# 
# outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_RP1_cells/count_mat_from_seurat.k36me3.rds"
# outf.stringent <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_RP1_cells/count_mat_from_seurat_stringent.k36me3.rds"
# 
# if (writeoutputs){
#   saveRDS(dat.sub.forLDA, file = outf)
#   saveRDS(dat.sub.forLDA.stringent, file = outf.stringent)
# }


# DimPlot(s.k562, reduction = "umap", cols = s.k562$col)

# 
# # Get indx ----------------------------------------------------------------
# 
# jmark.facs <- "K36"
# indir.facs <- file.path(hubpath, "jyeung/data/scChiC/facs_index_data/RP1_FUCCI_2020-09-25")
# infs.facs <- list.files(indir.facs, pattern = "*.csv", full.names = TRUE)
# print9(nfs.facs)
# 
# dats.clean <- lapply(infs, function(inf){
#   indx <- fread(inf)
#   bname <- strsplit(basename(inf), ".csv")[[1]][[1]]
#   
#   colnames(indx)
#   
#   grep("561", colnames(indx), value = TRUE)
#   cname.y <- "*[561] 585/29"
#   cname.x <- "*[488] 530/40"
#   
#   cnames.keep <- grepl(pattern = "585/29|530/40|Well", colnames(indx))
#   dat.clean <- indx[, ..cnames.keep]
#   dat.clean$bname <- bname
#   return(dat.clean)
# }) %>%
#   bind_rows()
# 
# ggplot(dats.clean, aes(y = `*[561] 585/29`, x = `*[488] 530/40`)) + 
#   geom_point(alpha = 0.2)  + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_x_log10() + scale_y_log10() 


