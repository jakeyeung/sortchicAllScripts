# Jake Yeung
# Date of Creation: 2021-01-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/4-make_cellcell_correlations.TSS_all.R
# 

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

# jmarks <- c("H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmark1 <- "H3K27me3"
jmark2 <- "H3K9me3"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K27me3_H3K9me3_analysis"
outpdf <- file.path(outdir, paste0("clustering_by_bins.", Sys.Date(), ".pdf"))
pdf(file = outpdf, useDingbats = FALSE)

# Load .metasmetas  -------------------------------------------------------------




ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned")
assertthat::assert_that(dir.exists(indir.meta))
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  print(inf)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})

dat.metas.reordered <- lapply(jmarks, function(jmark){
  if (jmark != "H3K9me3"){
    jdat.reordered <- dat.metas[[jmark]] 
    jdat.reordered$cluster <- factor(jdat.reordered$cluster, levels = ctypes)
    jdat.reordered <- jdat.reordered %>%
      arrange(cluster)
  } else {
    jdat.reordered <- dat.metas[[jmark]]
  }
  return(jdat.reordered)
})


cells.ordered.lst <- lapply(jmarks, function(jmark){
  dat.metas.reordered[[jmark]]$cell
})




# Load LDA ----------------------------------------------------------------



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

out.lda.lst <- lapply(outs.all.lst, function(jout) jout$out.lda)


count.mat.lst <- lapply(outs.all.lst, function(jout){
  jmat <- jout$count.mat
  jmat <- sweep(jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/")
  # jmat <- BinarizeMatrix(jmat)
})

tm.result.lst <- lapply(out.lda.lst, function(jout) AddTopicToTmResult(posterior(jout)))


dat.imputed.lst <- lapply(jmarks, function(jmark){
  jmat <- t(log2(tm.result.lst[[jmark]]$topic %*% tm.result.lst[[jmark]]$term))
  # reorder
  cells.ordered <- dat.metas.reordered[[jmark]]$cell
  jmat <- jmat[, cells.ordered]
})





# By topics ---------------------------------------------------------------


jmark <- "H3K27me3"

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

# labels_colors(dend) 

# jlabs <- labels(dend)
jlabs <- rownames(topics.mat)
jcols <- sapply(jlabs, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
colors_to_use <- jcols[order.dendrogram(dend)]

labels_colors(dend) <- colors_to_use
labels_colors(dend)

dend <- set(dend, "labels_cex", 0.5)

plot(dend, main = jmark)






# Load fits ---------------------------------------------------------------

jmark <- "H3K27me3"
inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
# inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.", jmark, ".2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData")

load(inf.fits, v=T)

params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
  mutate(log2fc = estimate / log(2))
params.long$padj <- p.adjust(params.long$pval.param)
jnames <- names(jfits.lst); names(jnames) <- jnames
pvals.long <- lapply(jnames, function(jname){
  x <- jfits.lst[[jname]]
  xvec <- x$pval
  data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

ggplot(params.long %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) + 
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Get DE genes  -----------------------------------------------------------


jmark.ref <- "H3K27me3"
# filter by pvals
# bins.filt <- subset(pvals.long.lst[[jmark.ref]], pval < 10^-10)$bin
bins.filt <- subset(pvals.long, pval < 10^-10)$bin

# Get imputed -------------------------------------------------------------


dat.imputed.lst <- lapply(jmarks, function(jmark){
  jmat <- t(log2(tm.result.lst[[jmark]]$topic %*% tm.result.lst[[jmark]]$term))
  # reorder
  cells.ordered <- dat.metas.reordered[[jmark]]$cell
  jmat <- jmat[, cells.ordered]
  # if (jmark == "H3K27me3"){
  #   rownames(jmat) <- sapply(rownames(jmat), function(rname) AssignHash(x = rname, jhash = bins.hash, null.fill = rname))
  # }
  return(jmat)
})




# Take top bins and make heatmap  -----------------------------------------

params.ref <- params.long

jparams <- unique(params.ref$param); names(jparams) <- jparams

# jparam <- "ClusterGranulocytes.Estimate"
params.annot.lst <- lapply(jparams, function(jparam){
  print(jparam)
  params.long.annot <- params.ref %>%
    filter(bin %in% bins.filt) %>%
    rowwise() %>%
    mutate(is.ctype = param == jparam) %>%
    group_by(bin, is.ctype) %>%
    filter(abs(log2fc) < 10) %>%
    # mutate(log2fc = ifelse(log2fc < -5, -5, log2fc),
    #        log2fc = ifelse(log2fc > 5, 5, log2fc)) %>%
    summarise(log2fc.mean = mean(log2fc)) %>%
    group_by(bin)  %>%
    filter(length(log2fc.mean) == 2) %>%
    summarise(log2fc.diff = log2fc.mean[[2]] - log2fc.mean[[1]]) %>%
    # arrange(desc(log2fc.diff)) %>%
    arrange(log2fc.diff) %>%
    mutate(param = jparam)
})

# take top
keepn <- 150

bins.top.filt.lst <- lapply(params.annot.lst, function(jdat){
  jdat <- jdat %>%
    # arrange(log2fc.diff)
  arrange(desc(log2fc.diff))
  jdat$bin[1:keepn]
})


# Make clustering ---------------------------------------------------------

jmark <- "H3K27me3"

cell2clst <- hash::hash(dat.metas[[jmark]]$cell, dat.metas[[jmark]]$cluster)
clst2col <- hash::hash(dat.metas[[jmark]]$cluster, dat.metas[[jmark]]$clustercol)

# topics.mat <- tm.result.lst[[jmark]]$topics
inmat <- t(dat.imputed.lst[[jmark]][unlist(bins.top.filt.lst), ])

rownames(inmat) <- sapply(rownames(inmat), function(x) cell2clst[[x]])

# rename 


jdist <- dist(inmat, method = "euclidean")
jclst <- hclust(jdist, method = "ward.D2")

dend <- as.dendrogram(jclst)

labels_colors(dend) 

# jlabs <- labels(dend)
jlabs <- rownames(inmat)
jcols <- sapply(jlabs, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
colors_to_use <- jcols[order.dendrogram(dend)]

labels_colors(dend) <- colors_to_use
labels_colors(dend)

dend <- set(dend, "labels_cex", 0.5)

plot(dend, main = jmark)


jlabs2 <- sapply(jlabs, function(x) Ctype2Lineage(x))
jcols2 <- sapply(jlabs2, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
bars_to_use <- jcols[order.dendrogram(dend)]
bars_to_use2 <- jcols2[order.dendrogram(dend)]



colored_bars(colors = bars_to_use, dend = dend, sort_by_labels_order = FALSE, rowLabels = "celltype") 
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

par(mar = c(10, 2, 1, 1))
dend %>% plot
plot(dend, main = jmark)
colored_bars(colors = bars_to_use2, dend = dend, sort_by_labels_order = FALSE, rowLabels = "celltype") 
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)



# shopw colored bars
par(mar = c(10, 2, 1, 1))
labels(dend) <- rep(".", nrow(inmat))
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


dev.off()



