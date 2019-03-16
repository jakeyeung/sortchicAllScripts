# Jake Yeung
# Date of Creation: 2019-03-06
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/canonical_correlation_bins.R
# On bins 

# https://satijalab.org/seurat/immune_alignment.html
# https://satijalab.org/seurat/Seurat_AlignmentTutorial.html

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

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load active marks LDA output --------------------------------------------



jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
jcolvec <- c("blue", "gray80", "red")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
names(out.objs) <- jmarks
names(out.objs.nobin) <- jmarks

tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))

mat.impute.lst <- lapply(tm.result.lst, function(tm.result) t(tm.result$topic %*% tm.result$term))

print(lapply(mat.impute.lst, dim))

# Wrangle data ------------------------------------------------------------

top.thres <- 0.999
top.regions.lst <- lapply(jmarks, function(jmark){
  topic.regions <- out.objs[[jmark]]$topic.regions  # for each cluster
  top.regions <- unique(unlist(topic.regions))
  return(top.regions) 
})
top.regions.merge <- unique(unlist(top.regions.lst))

mat.merged.lst <- lapply(mat.impute.lst, function(mat.impute){
  row.i <- which(top.regions.merge %in% rownames(mat.impute))
  return(as.data.frame(mat.impute[row.i, ]))
})

rnames.lst <- lapply(mat.merged.lst, function(x) rownames(x))
top.regions.intersect <- purrr::reduce(.x = rnames.lst, .f = intersect)

mat.merged.lst <- lapply(mat.merged.lst, function(mat.merged){
  return(mat.merged[top.regions.intersect, ])
})


# Do penalized PMD? -------------------------------------------------------


library(Seurat)

X <- as.matrix(mat.merged.lst[[1]])
Y <- as.matrix(mat.merged.lst[[2]])

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

hvg.k4me1 <-k4me1@assays$scChIC@var.features[1:1000]
hvg.k4me3 <-k4me3@assays$scChIC@var.features[1:1000]
hvg.union <- union(x = hvg.k4me1, y = hvg.k4me3)

marks.combined <- RunCCA(k4me1, k4me3, genes.use = hvg.union, num.cc = 30)

p1 <- DimPlot(object = marks.combined, group.by = "mark", reduction = "cca", pt.size = 0.5, dims = 3:4)
p2 <- VlnPlot(object = marks.combined, features = "CC_1", group.by = "mark")
plot_grid(p1, p2)



# print peak


jmark <- "H3K4me1"

# m1 <- PlotImputedPeaks2(tm.result.lst[[1]], jpeak, jmark,
#                         use.count.mat = NULL,
#                         usettings=dat.umap, 
#                         gname = jgene,
#                         jsize = jsize, jcolvec = jcolvec, .log = TRUE)

DimHeatmap(object = marks.combined, reduction = "cca", cells = 500, dim = 1:9, balanced = TRUE)
DimHeatmap(object = marks.combined, reduction = "cca", cells = 500, dim = 10:18, balanced = TRUE)



# marks.combined <- AlignSubspace(marks.combined, reduction.type = "cca", grouping.var = "mark", 
#                                  dims.align = 1:13)

# Plot top hits?   --------------------------------------------------------

nn.vec <- c(40, 35)
jmindist.vec <- c(0.2, 0.1)
jmetric <- "euclidean"
jseed <- 123
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed), 
                              nn.vec, jmindist.vec, SIMPLIFY = FALSE)
# plot umaps to check
topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks


# explore top peaks that are similar
print(head(sort(marks.combined@reductions$cca@feature.loadings[, 1], decreasing = TRUE)))
print(head(sort(marks.combined@reductions$cca@feature.loadings[, 2], decreasing = TRUE)))
print(head(sort(marks.combined@reductions$cca@feature.loadings[, 2], decreasing = FALSE)))

# jpeak <- "chr1:156920000-157020000"  # highest correlation
# jpeak <- "chr2:126420000-126520000"  # highest negative correlation
# jpeak <- head(names(sort(marks.combined@reductions$cca@feature.loadings[, 1], decreasing = TRUE)))[[5]]
jpeak <- head(names(sort(marks.combined@reductions$cca@feature.loadings[, 3], decreasing = FALSE)))[[1]]
jpeak <- head(names(sort(marks.combined@reductions$cca@feature.loadings[, 4], decreasing = TRUE)))[[1]]
jscale.fac <- 10^6
jpseudo <- 1
jgene <- "TopCC"
jmark <- "H3K4me1"

jcolvec <- c("blue", "gray80", "red")
m1 <- PlotImputedPeaks2(tm.result.lst[[1]], jpeak, jmark,
                        use.count.mat = NULL,
                        usettings=dat.umap.lst[[1]],
                        gname = jgene,
                        jsize = jsize, jcolvec = jcolvec, .log = TRUE)
m2 <- PlotImputedPeaks2(tm.result.lst[[2]], jpeak, jmark,
                        use.count.mat = NULL,
                        usettings=dat.umap.lst[[2]],
                        gname = jgene,
                        jsize = jsize, jcolvec = jcolvec, .log = TRUE)
multiplot(m1, m2, cols = 2)

PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = list(NULL, NULL), dat.umap.lst, 
                 jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)


# Plot the CCs ------------------------------------------------------------

plot(marks.combined@reductions$cca@feature.loadings[, 1], marks.combined@reductions$cca@feature.loadings[, 2], pch=21)
text(marks.combined@reductions$cca@feature.loadings[, 1], marks.combined@reductions$cca@feature.loadings[, 2], 
     labels = rownames(marks.combined@reductions$cca@feature.loadings))

# do UMAP 
cca.umap <- umap(marks.combined@reductions$cca@feature.loadings[, 1:10])
plot(cca.umap$layout[, 1], cca.umap$layout[, 2], pch = 20)
text(cca.umap$layout[, 1], cca.umap$layout[, 2], labels = rownames(cca.umap$layout))

# plot genes on far right
# jpeak <- names(sort(cca.umap$layout[, 1], decreasing = TRUE))[[2]]
jpeak <- names(sort(cca.umap$layout[, 2], decreasing = FALSE))[[1]]
# jpeak <- names(sort(cca.umap$layout[, 1]))[[2]]

m1 <- PlotImputedPeaks2(tm.result.lst[[1]], jpeak, jmark,
                        use.count.mat = NULL,
                        usettings=dat.umap.lst[[1]],
                        gname = jgene,
                        jsize = jsize, jcolvec = jcolvec, .log = TRUE)
m2 <- PlotImputedPeaks2(tm.result.lst[[2]], jpeak, jmark,
                        use.count.mat = NULL,
                        usettings=dat.umap.lst[[2]],
                        gname = jgene,
                        jsize = jsize, jcolvec = jcolvec, .log = TRUE)
multiplot(m1, m2, cols = 2)


