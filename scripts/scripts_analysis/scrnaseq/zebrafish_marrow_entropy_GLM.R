# Jake Yeung
# Date of Creation: 2019-08-16
# File: ~/projects/scchic/scripts/scripts_analysis/scrnaseq/zebrafish_marrow.R
# Zebrafish

#' ---
#' title: Re-analysis of Zebrafish kidney marrow scRNA-seq 
#' author: Jake Yeung
#' path: "~/projects/scchic/scripts/scripts_analysis/scrnaseq/zebrafish_marrow_entropy_GLM.R"
#' date: 2019-08-19
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---




#' ## Introduction
#' 
#' I am reanalyzing count tables sent by Chloe from Kidney Marrow. We want to overlay the "focusing" information onto the analysis. A natural measure of this focusing is entropy.
#' 
#' I am also taking this opportunity to try out better scRNA-seq analyses, namely avoiding the log transform and dividing the data by cell-specific library size. 
#' One way to avoid these procedures is to use a generalized linear model framework and model the mean-variance relationship using a negative binomial to remove the count noise.
#' Residuals from this noise model are assumed to be coming from biological variability. I think this results in a more robust analysis and finds more "structure" in the data compared to doing the log-transform and scalinig.
#' 
#' In this analysis, we find HSPCs to have higher entropy than other cell types, with the exception of thrombocytes (recall high entropy means more evenly distributed gene expression, low entropy means more focused on a few genes). 
#' This high entropy in thrombocytes might be a robust finding? A recent integrative analysis of mouse hematopoiesis found megakaryocytes (precursors of thrombocytes) have large number of cis-regulatory elements in the genome, comparable to HSCs (by counting peaks in ChIP-seq data), Xiang, Keller, biorxiv 2019.


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
library(here)

setwd(here())




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

GetGeneCounts <- function(x){
  # get number of zeros in x
  return(Matrix::nnzero(x) )
  
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
# mat.filt <- zf@assays$RNA@counts[sigVariedGene, ]



#' ## Variable gene selection (used to subselect for entropy calculation)

#' We still look for variable features for the entropy calculation, removes genes lowly expressed in all cells

mat <- zf@assays$RNA@counts
cells.keep <- which(colSums(mat) > 1000)
genes.keep <- which(rowMeans(mat) > 0)
mat <- mat[genes.keep, cells.keep]
gene.mean <- rowMeans(mat)
gene.var <- apply(mat, 1, function(jrow) var(jrow))
gene.cv2 <- gene.var / gene.mean ^ 2

# plot
print(range(gene.mean))
print(range(gene.cv2))

smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20, nrpoints = 0)
abline(a = 0, b = -1)


# plot(density(log10(Matrix::rowSums(mat))))
# plot(density(log10(Matrix::colSums(mat))))

minMeanForFit <- 10^-3
useForFit <- gene.mean >= minMeanForFit



# smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20)
# points(x = log10(gene.mean[useForFit]), y = log10(gene.cv2[useForFit]), pch = 20, col = 'red')

fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/gene.mean[useForFit] ), gene.cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients

xg <- 10^(seq( min(log10(gene.mean)), max(log10(gene.mean)), length.out=1000 ))
vfit <- a1/xg + a0

smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20, nrpoints = 0)
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

smoothScatter(x = log10(gene.mean), y = log10(gene.cv2), pch = 20, nrpoints = 0)
lines( log10(xg), log10(vfit), col="black", lwd=3 )
points(x = log10(gene.mean[sigVariedGenes]), y = log10(gene.cv2[sigVariedGenes]), pch = 20, col = 'red')


#' ## Calculate entropy

# Calculate entropy -------------------------------------------------------

mat.filt <- zf@assays$RNA@counts[sigVariedGenes, ]
S.vec <- apply(mat.filt, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
jmeta <- data.frame(S = S.vec)
rownames(jmeta) <- names(S.vec)

meta.merge <- left_join(meta, data.frame(rname = rownames(jmeta), jmeta)) %>%
  filter(!is.na(S))

m <- ggplot(meta.merge, aes(x = S, group = celltype)) + facet_wrap(~celltype, ncol = 1) + geom_density(fill = "lightblue") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Entropy [bits]")
print(m)

#' We find HSPCs tend to have higher entropy than other cell types. Exception is thrombocytes. 

#' ## Project entropy information onto UMAP
#' 
#' The input data matrix is the residuals of the technical noise model. Assuming deviations from the technical noise are coming from biological variability, then 
#' these residuals represent the variability of biological signal across the cell types. Hopefully they recapitulate Baron et al. and maybe find cell types that were not there before.
#' All cell type annotations were taken from Baron et al., but the UMAP comes from a different analysis. There might be a few cells defined as one cell type that wanders into a new cluster in this analysis, but the general picture should stil be correct.



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

#' If we are very keen we could even subdivide some of these annotated cell types into more fine-grained cell types? 
FeaturePlot(object = zf, features = 'entropy') + scale_color_viridis_c(end = 0.9) + ggtitle("Entropy [bits]")

