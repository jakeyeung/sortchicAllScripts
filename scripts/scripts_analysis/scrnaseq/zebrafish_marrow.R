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

inf <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.csv.gz"
inf.meta <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/tsne_clusterID_zebrafish_GateID_dataset.csv"

dat <- fread(inf)
meta <- fread(inf.meta)

colnames(meta) <- c("rname", "V1", "V2", "ClusterID", "Types")
colnames(dat)[[1]] <- "gene"

mat <- Matrix(as.matrix(as.data.frame(dat[, -1])), sparse = TRUE)
rownames(mat) <- dat$gene



# select variable features
gene.mean <- rowMeans(mat)
gene.var <- apply(mat, 1, function(jrow) var(jrow))
gene.cv2 <- gene.var / gene.mean ^ 2

# plot
smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20)
abline(a = 0, b = -1)


plot(density(log10(Matrix::rowSums(mat))))
plot(density(log10(Matrix::colSums(mat))))


# zf <- CreateSeuratObject(counts = mat[genes.remove, ])
# zf <- SCTransform(zf, verbose=TRUE)

# 
# 
# # fit curve
# minMeanForFit <- unname( quantile( gene.mean[ which( gene.cv2 > 0 ) ], .95 ) )
# minMeatForFit <- 10^-2
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

plot(density(log10(Matrix::rowSums(mat))))
plot(density(log10(Matrix::colSums(mat))))

cells.keep <- which(colSums(mat) > 1000)

plot(density(log10(colSums(mat[, cells.keep]))))

# try Seurat


zf <- CreateSeuratObject(counts = mat, project = "zf", min.cells = 3, min.features = 200)

# plot1 <- FeatureScatter(zf, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)

zf <- FindVariableFeatures(zf, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(zf), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(zf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# add labels?
zf@meta.data$cell <- rownames(zf@meta.data)
jhash <- hash(meta$rname, meta$Types)
jhash2 <- hash(meta$rname, meta$ClusterID)
zf@meta.data$celltype <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash[[x]], NA)
})
zf@meta.data$clusterid <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash2[[x]]), jhash2[[x]], NA)
})

all.genes <- rownames(zf)
zf <- ScaleData(zf, features = all.genes)
zf <- RunPCA(zf, features = VariableFeatures(object = zf))
zf <- RunUMAP(zf, dims = 1:10)
DimPlot(zf, reduction = "umap", label = FALSE, group.by = "clusterid", na.value = "#C0C0C0") + scale_color_viridis_d()

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
