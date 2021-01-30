# Jake Yeung
# Date of Creation: 2021-01-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/4-make_cellcell_correlations.R
# Explore hierarchical structure in K27me3 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(topicmodels)

library(dendextend)


# Functions ---------------------------------------------------------------


# ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")

Ctype2Lineage <- function(x){
  if (x == "NKs"){
    xnew <- "Bcells"
  } else if (x == "Basophils"){
    xnew <- "Granulocytes"
  } else if (x == "pDCs"){
    xnew <- "Bcells"
  } else if (x == "DCs"){
    xnew <- "Granulocytes"
  } else {
    xnew <- x
  }
  return(xnew)
}

hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load LDAs ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# load LDA 
outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

tm.result.lst <- lapply(jmarks, function(jmark){
  tm.result <- posterior(outs.all.lst[[jmark]]$out.lda)
  tm.result <- AddTopicToTmResult(tm.result, jsep = "")
  return(tm.result)
})


# Load .metasmetas  -------------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})



# Do cellcell distances?  -------------------------------------------------



# jmark <- "H3K4me3"

jmark <- "H3K9me3"
jmark <- "H3K27me3"


jmark <- "H3K4me1"

jmark <- "H3K4me3"

cell2clst <- hash::hash(dat.metas[[jmark]]$cell, dat.metas[[jmark]]$cluster)
clst2col <- hash::hash(dat.metas[[jmark]]$cluster, dat.metas[[jmark]]$clustercol)



topics.mat <- tm.result.lst[[jmark]]$topics
# transform
# topics.mat <- log2(topics.mat / (1 - topics.mat))

rownames(topics.mat) <- sapply(rownames(topics.mat), function(x) cell2clst[[x]])

# rename 


jdist <- dist(topics.mat, method = "euclidean")
jclst <- hclust(jdist, method = "ward.D2")

dend <- as.dendrogram(jclst)

labels_colors(dend) 

# jlabs <- labels(dend)
jlabs <- rownames(topics.mat)
jcols <- sapply(jlabs, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
colors_to_use <- jcols[order.dendrogram(dend)]

labels_colors(dend) <- colors_to_use
labels_colors(dend)

dend <- set(dend, "labels_cex", 0.5)

plot(dend, main = jmark)



# Try clustering adjusted impuetd mat -------------------------------------

# jmark <- "H3K27me3"
# jmark <- "H3K4me3"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K27me3_H3K9me3_analysis/clustering_outputs_marks.", Sys.Date(), ".lineage.pdf")
pdf(outpdf, useDingbats = FALSE)

for (jmark in jmarks){
  print(jmark)
  
  if (jmark != "H3K9me3"){
    inf.mat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.", jmark, ".2021-01-08.rearranged.RData")
    load(inf.mat, v=T)  # mat.adj.tmp
    jmat.input <- t(mat.adj.tmp)
  } else {
    inf.mat2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2020-12-20.H3K27me3rep2rep3reseq.RData"
    load(inf.mat2, v=T)
    jmat.input.tmp <- mat.adj.lst$H3K9me3
    rnames.tmp <- mat.adj.lst$H3K9me3$rname
    jmat.input.tmp$rname <- NULL
    jmat.input.tmp <- as.matrix(jmat.input.tmp)
    rownames(jmat.input.tmp) <- rnames.tmp
    jmat.input <- t(jmat.input.tmp)
  }
  
  cell2clst <- hash::hash(dat.metas[[jmark]]$cell, dat.metas[[jmark]]$cluster)
  clst2col <- hash::hash(dat.metas[[jmark]]$cluster, dat.metas[[jmark]]$clustercol)
  
  
  rownames(jmat.input) <- sapply(rownames(jmat.input), function(x) AssignHash(x = x, jhash = cell2clst, null.fill = x))
  
  jdist <- dist(jmat.input, method = "euclidean")
  jclst <- hclust(jdist, method = "ward.D2")
  
  dend <- as.dendrogram(jclst)
  
  labels_colors(dend) 
  
  jlabs <- rownames(jmat.input)
  jlabs2 <- sapply(jlabs, function(x) Ctype2Lineage(x))
  jcols <- sapply(jlabs, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
  jcols2 <- sapply(jlabs2, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
  colors_to_use <- jcols[order.dendrogram(dend)]
  bars_to_use <- jcols[order.dendrogram(dend)]
  bars_to_use2 <- jcols2[order.dendrogram(dend)]
  
  labels_colors(dend) <- colors_to_use
  labels_colors(dend)
  
  # dend <- set(dend, "labels_cex", 0.5)
  # plot(dend, main = jmark)
  
  # shopw colored bars
  par(mar = c(10, 2, 1, 1))
  labels(dend) <- rep(".", nrow(jmat.input))
  # dend %>% set("labels_to_char") %>% set("labels", rep(".", nrow(jmat.input))) %>% plot
  dend %>% plot
  plot(dend, main = jmark)
  colored_bars(colors = bars_to_use, dend = dend, sort_by_labels_order = FALSE, rowLabels = "celltype") 
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  
  par(mar = c(10, 2, 1, 1))
  dend %>% plot
  plot(dend, main = jmark)
  colored_bars(colors = bars_to_use2, dend = dend, sort_by_labels_order = FALSE, rowLabels = "celltype") 
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  
}
dev.off()


