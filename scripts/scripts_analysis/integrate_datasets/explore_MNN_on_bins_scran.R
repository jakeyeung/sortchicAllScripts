# Jake Yeung
# Date of Creation: 2019-02-22
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/explore_MNN_on_bins_scran.R
# Use scran package

# rm(list=ls())

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

# merge all
# mat.merged <- purrr::reduce(.x = mat.merged.lst, .f = dplyr::bind_cols)  # silently drop unmatched names?? No then manually intersectjA
# rownames(mat.merged) <- top.regions.intersect

# merge active marks only
# mat.merged <- purrr::reduce(.x = mat.merged.lst, .f = dplyr::bind_cols)  # silently drop unmatched names?? No then manually intersect
mat.merged <- bind_cols(mat.merged.lst[["H3K4me1"]], mat.merged.lst[["H3K4me3"]])
rownames(mat.merged) <- top.regions.intersect

coldat <- data.frame(tech = sapply(colnames(mat.merged), function(x) strsplit(x, "_")[[1]][[2]]))
rownames(coldat) <- colnames(mat.merged)

bm <- CreateSeuratObject(counts = mat.merged, meta.data = coldat)
bm.lst <- SplitObject(object = bm, split.by = "tech")

# get variable genes
for (i in 1:length(x = bm.lst)) {
  bm.lst[[i]] <- FindVariableFeatures(object = bm.lst[[i]], selection.method = "dispersion", nfeatures = 500, verbose = TRUE)
}

# sort(jmarks.all)
reference.list <- bm.lst[sort(unique(coldat$tech))]

bm.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

bm.integrated <- IntegrateData(anchorset = bm.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = bm.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
bm.integrated <- ScaleData(object = bm.integrated, verbose = FALSE)
bm.integrated <- RunPCA(object = bm.integrated, npcs = 100, verbose = FALSE)
bm.integrated <- RunUMAP(object = bm.integrated, reduction = "pca", 
                               dims = 1:100)
p1 <- DimPlot(object = bm.integrated, reduction = "umap", group.by = "tech")
print(p1)

