# Jake Yeung
# Date of Creation: 2019-03-06
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/canonical_correlation.R
# Canonical correlation to integrate datasets



rm(list=ls())

setwd("~/projects/scchic")

library(GGally)
library(purrr)

library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

library(cowplot)



library(GGally)

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")




# Get dirs ----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3")
mdirs <- lapply(jmarks, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})
jmarks.repress <- c("H3K27me3", "H3K9me3")
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  print(mdir)
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})
jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all
mdirs <- c(mdirs, mdirs.repress)
names(mdirs) <- c(jmarks, jmarks.repress)
mara.outs <- lapply(mdirs, LoadMARA)
head(mara.outs[[1]]$act.long)


# CCA analysis ------------------------------------------------------------

# do on activation marks
X <- as.matrix(mara.outs$H3K4me1$act.mat %>% dplyr::select(-motif)); rownames(X) <- mara.outs$H3K4me1$act.mat$motif
Y <- as.matrix(mara.outs$H3K4me3$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs$H3K4me3$act.mat$motif
# Z <- cancor(X, Y, xcenter = TRUE, ycenter = TRUE)

# need to do penalized because X and Y are low rank
# source("~/projects/rmscca/scca_CVperm.R")
# source("~/projects/rmscca/scca_function.R")
# source("~/projects/rmscca/sample_sigma12_function.R")
library(PMA)

dat <- list(X = X, Y = Y)
# out <- PMA::MultiCCA(dat, ncomponents = 2)
out <- PMA::CCA(x = X, z = Y, K = 2)

plot(out$u[, 1], out$u[, 2])
plot(out$v[, 1], out$v[, 2])

X.cca <- X %*% out$u
Y.cca <- Y %*% out$v

par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(X.cca[, 1], X.cca[, 2])
text(X.cca[, 1], X.cca[, 2], labels = mara.outs$H3K4me1$act.mat$motif)
plot(Y.cca[, 1], Y.cca[, 2])
text(Y.cca[, 1], Y.cca[, 2], labels = mara.outs$H3K4me3$act.mat$motif)

# which cells show correlation?

# plot cells with non-zero weights onto the UMAP



# plot(out$ws[[1]][, 1], out$ws[[1]][, 2])
# plot(out$ws[[2]][, 1], out$ws[[2]][, 2])


# Use Seurat implementation -----------------------------------------------

# https://satijalab.org/seurat/immune_alignment.html


# Seurat 3.0 implementation -----------------------------------------------


library(Seurat)

jmeta1 <- data.frame(cell = colnames(X), mark = "H3Km4e1", stringsAsFactors = FALSE)
jmeta3 <- data.frame(cell = colnames(Y), mark = "H3Km4e3", stringsAsFactors = FALSE)
k4me1 <- CreateSeuratObject(counts = X, meta.data = jmeta1, assay = "scChIC")
k4me3 <- CreateSeuratObject(counts = Y, meta.data = jmeta3, assay = "scChIC")

# k4me1 <- NormalizeData(k4me1, normalization.method = "RC")
k4me1 <- ScaleData(k4me1, do.scale = FALSE)
# k4me3 <- NormalizeData(k4me3, normalization.method = "RC")
k4me3 <- ScaleData(k4me3, do.scale = FALSE)

k4me1 <- FindVariableFeatures(k4me1, nfeatures = 1000)
k4me3 <- FindVariableFeatures(k4me3, nfeatures = 1000)

k4me1@meta.data$mark <- "H3K4me1"
k4me3@meta.data$mark <- "H3K4me3"

genes.use <- intersect(mara.outs$H3K4me3$act.mat$motif, mara.outs$H3K4me1$act.mat$motif)

marks.combined <- RunCCA(k4me1, k4me3, genes.use = genes.use, num.cc = 30)

p1 <- DimPlot(object = marks.combined, group.by = "mark", reduction = "cca", pt.size = 0.5, dims = 3:4)
p2 <- VlnPlot(object = marks.combined, features = "CC_1", group.by = "mark")
plot_grid(p1, p2)

# # PrintDim(object = marks.combined, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
# p3 <- MetageneBicorPlot(marks.combined, grouping.var = "mark", dims.eval = 1:30, 
#                         display.progress = FALSE)


# Seurat 2.4 implementation -----------------------------------------------

library(devtools)
dev_mode(T)

# devtools::install_cran("Seurat")
devtools::install_github(repo = "mojaveazure/loomR", ref = "master")
devtools::install_cran("diffusionMap")
devtools::install_cran("FNN")
devtools::install_cran("caret")
devtools::install_cran("tclust")
devtools::install_cran("ranger")
devtools::install_github(repo = 'satijalab/seurat@f4da9c6', dependencies = FALSE)  # v2.2.1

library(Seurat)

jmeta1 <- data.frame(cell = colnames(X), mark = "H3Km4e1", stringsAsFactors = FALSE)
jmeta3 <- data.frame(cell = colnames(Y), mark = "H3Km4e3", stringsAsFactors = FALSE)
k4me1 <- CreateSeuratObject(counts = X, meta.data = jmeta1, assay = "scChIC")
k4me3 <- CreateSeuratObject(counts = Y, meta.data = jmeta3, assay = "scChIC")

# k4me1 <- NormalizeData(k4me1, normalization.method = "RC")
k4me1 <- ScaleData(k4me1, do.scale = FALSE)
# k4me3 <- NormalizeData(k4me3, normalization.method = "RC")
k4me3 <- ScaleData(k4me3, do.scale = FALSE)

k4me1 <- FindVariableFeatures(k4me1, nfeatures = 1000)
k4me3 <- FindVariableFeatures(k4me3, nfeatures = 1000)

k4me1@meta.data$mark <- "H3K4me1"
k4me3@meta.data$mark <- "H3K4me3"

genes.use <- intersect(mara.outs$H3K4me3$act.mat$motif, mara.outs$H3K4me1$act.mat$motif)

marks.combined <- RunCCA(k4me1, k4me3, genes.use = genes.use, num.cc = 30)

p1 <- DimPlot(object = k4me1, reduction.use = "cca", group.by = "mark", pt.size = 0.5, 
              do.return = TRUE)
p2 <- VlnPlot(object = k4me3, features.plot = "CC1", group.by = "mark", do.return = TRUE)
plot_grid(p1, p2)

p1 <- DimPlot(object = marks.combined, group.by = "mark", reduction = "cca", pt.size = 0.5)
p2 <- VlnPlot(object = marks.combined, features = "CC_1", group.by = "mark")
plot_grid(p1, p2)

dev_mode(F)

