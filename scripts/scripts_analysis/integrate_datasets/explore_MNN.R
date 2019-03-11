# Jake Yeung
# Date of Creation: 2019-02-16
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/explore_MNN.R
# Integrate H3K4me1 and H3K4me3 using mutual nearest neighbors
# or look at the FACS data


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

# Functions ---------------------------------------------------------------




# Get dirs ----------------------------------------------------------------

plotf <- "~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/mara/all_4_marks_motifs.pdf"
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
mara.outs <- lapply(mdirs, LoadMARA)
head(mara.outs[[1]]$act.long)

# act.long.merged <- rbind(mara.outs[[1]]$act.long %>% mutate(mark = jmarks.all[[1]]),
#                          mara.outs[[2]]$act.long %>% mutate(mark = jmarks.all[[2]]), 
#                          mara.outs[[3]]$act.long %>% mutate(mark = jmarks.all[[3]]), 
#                          mara.outs[[4]]$act.long %>% mutate(mark = jmarks.all[[4]]))
# integrate the datasets together
# ggplot(act.long.merged, aes(x = activity)) + geom_density() + facet_wrap(~mark) + 
# theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long, mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)

# Run for only 2 of these 
# act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long)
act.long.merged <- rbind(mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)
# act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[3]]$act.long)
# act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[4]]$act.long)


# # check what data structure I need
# pancreas.data <- readRDS(file = "~/Downloads/pancreas_v3_files/pancreas_expression_matrix.rds")  # genes by cells
# metadata <- readRDS(file = "~/Downloads/pancreas_v3_files/pancreas_metadata.rds")  # rownames are cells, columns are "tech" and "celltype"
# pancreas <- CreateSeuratObject(counts = pancreas.data, meta.data = metadata)
# pancreas.list <- SplitObject(object = pancreas, split.by = "tech")
# for (i in 1:length(x = pancreas.list)) {
#   pancreas.list[[i]] <- NormalizeData(object = pancreas.list[[i]], verbose = FALSE)
#   pancreas.list[[i]] <- FindVariableFeatures(object = pancreas.list[[i]], 
#                                              selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# }
# 
# 

act.mat <- spread(act.long.merged, key = cell, value = activity)
rownames(act.mat) <- act.mat$motif; act.mat$motif <- NULL

coldat <- data.frame(tech = sapply(colnames(act.mat), function(x) strsplit(x, "_")[[1]][[2]]))
rownames(coldat) <- colnames(act.mat)

bm <- CreateSeuratObject(counts = act.mat, meta.data = coldat)
bm.lst <- SplitObject(object = bm, split.by = "tech")
# get variable genes
for (i in 1:length(x = bm.lst)) {
  bm.lst[[i]] <- FindVariableFeatures(object = bm.lst[[i]], selection.method = "dispersion", nfeatures = 300, verbose = TRUE)
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
                         dims = 1:30)
p1 <- DimPlot(object = bm.integrated, reduction = "umap", group.by = "tech")
print(p1)


# Go deeper: what is the anchor cells and stuff/ --------------------------


