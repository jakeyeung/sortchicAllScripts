# Jake Yeung
# Date of Creation: 2019-08-16
# File: ~/projects/scchic/scripts/scripts_analysis/scrnaseq/zebrafish_marrow.R
# Zebrafish



rm(list=ls())

require(statmod)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(irlba)
library(umap)
library(scchicFuncs)
library(Seurat)
library(hash)
library(xlsx)

# Functions ---------------------------------------------------------------

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


# Load  -------------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.csv.gz"
inf.meta <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/tsne_clusterID_zebrafish_GateID_dataset.csv"
inf.meta2 <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/cluster_info_JY_edited.xlsx"

dat <- fread(inf, stringsAsFactors = FALSE)
meta <- fread(inf.meta, stringsAsFactors = TRUE)
meta2 <- read.xlsx2(inf.meta2, sheetIndex = 1, header = FALSE, colClasses = c("integer", "character")); colnames(meta2) <- c("ClusterID", "celltype")


colnames(meta) <- c("rname", "V1", "V2", "ClusterID", "experi")
colnames(dat)[[1]] <- "gene"

meta <- left_join(meta, meta2)

# add celltype from excle fie 





# Load GLM output ---------------------------------------------------------


# run a negbinom regression to estimate expression and remove count noise (analogous to DESeq2)
zf <- readRDS("/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset_PenalizedNegBinRegress.rds")

# Calculate entropy -------------------------------------------------------

# # do it on the raw counts... also try the denoised counts
# genes.keep <- gsub("-", "_", rownames(zf@assays$RNA))  # make rownames compatible with raw mat
# cells.keep <- colnames(zf@assays$RNA)
# 
# genes.keep.i <- which(rownames(mat) %in% genes.keep)
# cells.keep.i <- which(colnames(mat) %in% cells.keep)

# mat.filt <- mat[sigVariedGene, cells.keep]
mat.filt <- zf@assays$RNA@counts[sigVariedGene, ]
# mat.filt <- sweep(exp(mat.filt), MARGIN = 2, STATS = colSums(exp(mat.filt)), FUN = "/")

S.vec <- apply(mat.filt, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
jmeta <- data.frame(S = S.vec)
rownames(jmeta) <- names(S.vec)

# Do umap -----------------------------------------------------------------



# add meta data
zf@meta.data$cell <- rownames(zf@meta.data)
jhash <- hash(meta$rname, meta$experi)
jhash2 <- hash(meta$rname, meta$ClusterID)
jhash3 <- hash(meta$rname, as.character(meta$celltype))
jhash.entropy <- hash(rownames(jmeta), jmeta$S)
zf@meta.data$experi <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash[[x]], NA)
})
zf@meta.data$clusterid <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash2[[x]]), jhash2[[x]], NA)
})
zf@meta.data$celltype <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash3[[x]]), jhash3[[x]], NA)
})
zf@meta.data$entropy <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash.entropy[[x]]), jhash.entropy[[x]], NA)
})

zf <- RunPCA(zf, verbose = FALSE)
zf <- RunUMAP(zf, dims = 1:30, verbose = FALSE)

zf <- FindNeighbors(zf, dims = 1:30, verbose = FALSE)
zf <- FindClusters(zf, verbose = FALSE)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
DimPlot(zf, label = TRUE, group.by = "celltype", cols = cbPalette)
FeaturePlot(object = zf, features = 'entropy') + scale_color_viridis_c()

# add the cell labels


# try my own variable gene implementation

# select variable features

mat <- mat.filt
gene.mean <- rowMeans(mat)
gene.var <- apply(mat, 1, function(jrow) var(jrow))
gene.cv2 <- gene.var / gene.mean ^ 2

# plot
smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20)
abline(a = 0, b = -1)


plot(density(log10(Matrix::rowSums(mat))))
plot(density(log10(Matrix::colSums(mat))))

minMeanForFit <- 10^-3
useForFit <- gene.mean >= minMeanForFit # & spikeins

print(length(which(useForFit)))

smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20)
points(x = log10(gene.mean[useForFit]), y = log10(gene.cv2[useForFit]), pch = 20, col = 'red')

fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/gene.mean[useForFit] ), gene.cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

xg <- 10^(seq( min(log10(gene.mean)), max(log10(gene.mean)), length.out=1000 ))
vfit <- a1/xg + a0

smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20)
lines( log10(xg), log10(vfit), col="black", lwd=3 )
# add CI
# dof <- nrow(mat) - 1
# dof <- 
dof <- ncol(mat) - 1
lines(log10(xg),log10(vfit * qchisq(0.95,dof)/dof),lty=1,col="black")
lines(log10(xg),log10(vfit * qchisq(0.05,dof)/dof),lty=1,col="black")

afit <- a1/gene.mean+a0
varFitRatio <- gene.var/(afit*gene.mean^2)

pval <- pchisq(varFitRatio*dof,df=dof,lower.tail=F)
adj.pval <- p.adjust(pval,"fdr")
sigVariedGenes <- adj.pval<1e-3;
table(sigVariedGenes)





# 
# 
# # Calculate entropy -------------------------------------------------------
# 
# 
# 
# 
# zf <- CreateSeuratObject(counts = mat[, cells.keep], project = "zf", min.cells = 3, min.features = 200, meta.data = jmeta)
# # zf <- CreateSeuratObject(counts = mat, project = "zf", min.cells = 3, min.features = 200, meta.data = jmeta)
# 
# # zf@assays$RNA@var.features <- names(which(sigVariedGenes))
# VariableFeatures(zf) <- sigVariedGenes
# 
# # zf <- FindVariableFeatures(zf, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(zf), 10)
# 
# # plot variable features with and without labels
# # plot1 <- VariableFeaturePlot(zf)
# # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# # CombinePlots(plots = list(plot1, plot2))
# 
# # add labels?
# zf@meta.data$cell <- rownames(zf@meta.data)
# jhash <- hash(meta$rname, meta$experi)
# jhash2 <- hash(meta$rname, meta$ClusterID)
# jhash3 <- hash(meta$rname, as.character(meta$celltype))
# zf@meta.data$experi <- sapply(zf@meta.data$cell, function(x){
#   ifelse(!is.null(jhash[[x]]), jhash[[x]], NA)
# })
# zf@meta.data$clusterid <- sapply(zf@meta.data$cell, function(x){
#   ifelse(!is.null(jhash2[[x]]), jhash2[[x]], NA)
# })
# zf@meta.data$celltype <- sapply(zf@meta.data$cell, function(x){
#   ifelse(!is.null(jhash3[[x]]), jhash3[[x]], NA)
# })
# 
# all.genes <- rownames(zf)
# zf <- ScaleData(zf, features = all.genes)
# zf <- RunPCA(zf, features = VariableFeatures(object = zf))
# zf <- RunUMAP(zf, dims = 1:10)
# DimPlot(zf, reduction = "umap", label = FALSE, group.by = "clusterid", na.value = "#C0C0C0") + scale_color_viridis_d()
# 
# # label by entropy
# # DimPlot(zf, reduction = "umap", label = FALSE, group.by = "clusterid", na.value = "#C0C0C0") + scale_color_viridis_d()
# 
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
# DimPlot(zf, reduction = "umap", label = FALSE, group.by = "celltype", na.value = "#C0C0C0", cols = cbPalette) 
# 
# FeaturePlot(object = zf, features = 'S') + scale_color_viridis_c()



# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")


# lsi.out <- scchicFuncs::RunLSI(t(as.matrix(mat[sigVariedGenes, cells.keep])))
#
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.15
# umap.out <- umap(as.matrix(lsi.out$v), config = jsettings)
# dat.umap.out <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
#
# ggplot(dat.umap.out, aes(x = umap1, y = umap2)) + geom_point() +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# add labels


# zf <- CreateSeuratObject(counts = mat[, cells.keep])
# zf <- SCTransform(zf, verbose=TRUE)

# pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
