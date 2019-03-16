# Jake Yeung
# Date of Creation: 2019-03-07
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/canonical_correlation_rewrite.R
# Redo canonical correlation analysis and allow opposite signs to happen downstream


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

jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"
# do on activation marks
X <- as.matrix(mara.outs[[jmark1]]$act.mat %>% dplyr::select(-motif)); rownames(X) <- mara.outs[[jmark1]]$act.mat$motif
# Y <- as.matrix(mara.outs$H3K4me3$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs$H3K4me3$act.mat$motif
Y <- as.matrix(mara.outs[[jmark2]]$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs[[jmark2]]$act.mat$motif
# Z <- cancor(X, Y, xcenter = TRUE, ycenter = TRUE)

# need to do penalized because X and Y are low rank
# source("~/projects/rmscca/scca_CVperm.R")
# source("~/projects/rmscca/scca_function.R")
# source("~/projects/rmscca/sample_sigma12_function.R")
# library(PMA)
# 
# dat <- list(X = X, Y = Y)
# # out <- PMA::MultiCCA(dat, ncomponents = 2)
# out <- PMA::CCA(x = X, z = Y, K = 2)
# 
# plot(out$u[, 1], out$u[, 2])
# plot(out$v[, 1], out$v[, 2])
# 
# X.cca <- X %*% out$u
# Y.cca <- Y %*% out$v
# 
# par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# plot(X.cca[, 1], X.cca[, 2])
# text(X.cca[, 1], X.cca[, 2], labels = mara.outs[[jmark1]]$act.mat$motif)
# plot(Y.cca[, 1], Y.cca[, 2])
# text(Y.cca[, 1], Y.cca[, 2], labels = mara.outs[[jmark2]]$act.mat$motif)

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

p1 <- DimPlot(object = marks.combined, group.by = "mark", reduction = "cca", pt.size = 0.5, dims = 2:3)
p2 <- VlnPlot(object = marks.combined, features = "CC_1", group.by = "mark")
plot_grid(p1, p2)


# Find correlated genes ---------------------------------------------------

loadings <- marks.combined@reductions$cca@feature.loadings
embeddings <- marks.combined@reductions$cca@cell.embeddings

# reproduce cell embeddings: should be diagonalized SVD of CCA

# Run the diagonal canonical correlation procedure
#
# @param mat1         First matrix
# @param mat2         Second matrix
# @param standardize  Standardize matrices - scales columns to have unit
#                     variance and mean 0
# @param k            Number of canonical correlation vectors (CCs) to calculate
#
# @return             Returns the canonical correlation vectors - corresponding
#                     to the left and right singular vectors after SVD - as well
#                     as the singular values.
#
library(irlba)
jCanonCor <- function(mat1, mat2, k = 20) {
  set.seed(seed = 42)
  # if (standardize) {
  #   mat1 <- Standardize(mat = mat1, display_progress = FALSE)
  #   mat2 <- Standardize(mat = mat2, display_progress = FALSE)
  # }
  # mat3 <- FastMatMult(m1 = t(x = mat1), m2 = mat2)
  mat3 <- t(mat1) %*% mat2
  cca.svd <- irlba::irlba(A = mat3, nv = k)
  return(list(u = cca.svd$u, v = cca.svd$v, d = cca.svd$d))
}

X <- k4me1@assays$scChIC@scale.data
Y <- k4me3@assays$scChIC@scale.data
cca.results <- jCanonCor(X, Y, k = 20)

cca.data <- rbind(cca.results$u, cca.results$v)

colnames(x = cca.data) <- paste0("CC", 1:20)
rownames(cca.data) <- c(colnames(X), colnames(Y))

# why does this matter? Maybe it does? 
cca.data.flip <- apply(cca.data, MARGIN = 2, function(x){
  if(sign(x[1]) == -1) {
    x <- x * -1
  }
  return(x)
})

# get minimum absolute number, return with actual sign
Vectorize(SelectAbsMin <- function(x1, x2){
  return(c(x1, x2)[[which.min(c(x1, x2))]])
}, vectorize.args = c("x1", "x2"), SIMPLIFY = FALSE, USE.NAMES = TRUE)


# how are these feature loadings calculated? Just projection
embeds <- cca.data.flip
# average across embeds
loads <- t(t(embeds) %*% rbind(t(X), t(Y)))
# take min function?
embeds1 <- cca.results$u
embeds2 <- cca.results$v

loads1 <- t(t(embeds1) %*% t(X))
loads2 <- t(t(embeds2) %*% t(Y))

loads.min <- matrix(mapply(SelectAbsMin, loads1, loads2), nrow = nrow(loads1), ncol = ncol(loads1), dimnames = list(rownames(loads1), colnames(loads1)))


# compare embeds and loads
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(marks.combined@reductions$cca@feature.loadings[, 1], marks.combined@reductions$cca@feature.loadings[, 2], pch = 20)
text(marks.combined@reductions$cca@feature.loadings[, 1], marks.combined@reductions$cca@feature.loadings[, 2], 
     labels = rownames(marks.combined@reductions$cca@feature.loadings))

plot(loads[, 1], loads[, 2], pch = 20)
text(loads[, 1], loads[, 2], 
     labels = rownames(loads))

plot(loads.min[, 1], loads.min[, 2], pch = 20)
text(loads.min[, 1], loads.min[, 2], 
     labels = rownames(loads.min))


# Plot hits ---------------------------------------------------------------

nn.vec <- c(40, 35)
jmindist.vec <- c(0.2, 0.1)
jmetric <- "euclidean"
jmindist=0.1
jseed=123
custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])

out.lda.lst <- lapply(infs[c(jmark1, jmark2)], LoadLDA)

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)
names(umap.lda.lst) <- c(jmark1, jmark2)
names(mara.outs) <- jmarks.all
dat.merged.lst <- lapply(c(jmark1, jmark2), function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- c(jmark1, jmark2)

jmotif <- "Hoxc6"

jmotif <- "Cebpb"

jmotif <- "Gata3"
jmotif <- "Gfi1"
jmotif <- "Hoxa2"
jmotif <- "Bcl3"

jmotif <- "Tal1"
jmotif <- "Zbtb16"
jmotif <- "Hoxa2"
jmotif <- "Hoxa1"
jmotif <- "Zfp110"
jmotif <- "Gfi1"
jmotif <- "Foxc1"
jmotif <- "Mecp2"
jmotif <- "Hic1"
jmotif <- "Atf2"
m1 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark1]], mara.outs[[jmark1]]$zscores, jmark1, jsize = 0.75)
m2 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark2]], mara.outs[[jmark2]]$zscores, jmark2, jsize = 0.75)
multiplot(m1, m2)

plot(loads.min[, 1], loads.min[, 2], pch = 20)
text(loads.min[, 1], loads.min[, 2],
     labels = rownames(loads.min))

# loads.umap <- umap(loads.min[, 1:10])
# plot(loads.umap$layout[, 1], loads.umap$layout[, 2], pch = 20)
# text(loads.umap$layout[, 1], loads.umap$layout[, 2], labels = rownames(loads.umap$layout))

# Now let’s do multiCCA ---------------------------------------------------



