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
outpdf <- file.path(outdir, paste0("clustering_by_TSS.", Sys.Date(), ".pdf"))
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

# jmark <- "H3K27me3"
# jmark2 <- "H3K9me3"

indir.lda.others <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows"
indir.lda.k27me3 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-H3K27me3-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE"

outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    fname.lda <- paste0("lda_outputs.count_mat_from_TSS.", jmark, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_10000.K-30.Robj")
    # inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
    inf.lda <- file.path(indir.lda.others, fname.lda)
  } else {
    # fname.lda <- paste0("lda_outputs.count_mat_from_TSS.", jmark, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_10000.K-30.Robj")
    fname.lda <- paste0("ldaOut.PZ-BM-rep3-", jmark, "-rep2rep3reseq.TSS.varfilt.K-30.Robj")
    inf.lda <- file.path(indir.lda.k27me3, fname.lda)
  }
  print(inf.lda)
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


# order topics: K27me3 only 
topics.ordered <- OrderTopicsByEntropy(tm.result.lst$H3K27me3)




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

labels_colors(dend) 

# jlabs <- labels(dend)
jlabs <- rownames(topics.mat)
jcols <- sapply(jlabs, function(x) AssignHash(x, jhash = clst2col, null.fill = "grey85"))
colors_to_use <- jcols[order.dendrogram(dend)]

labels_colors(dend) <- colors_to_use
labels_colors(dend)

dend <- set(dend, "labels_cex", 0.5)

plot(dend, main = jmark)




# Get DE from TSS ---------------------------------------------------------

outs.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.fits.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.", jmark, ".2020-12-12.newannot2.witherrors.TSS.RData"))
  
  load(inf.fits.tmp, v=T)
  
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
  
  return(list(params.long = params.long, pvals.long = pvals.long))
})

params.long.lst <- lapply(outs.lst, function(jout) jout$params.long)
pvals.long.lst <- lapply(outs.lst, function(jout) jout$pvals.long)






# Relabel H3K27me3  -------------------------------------------------------

library(hash)

bins.k27me3 <- params.long.lst$H3K27me3$bin
bins.others <- params.long.lst$H3K4me1$bin
bins.others.coords <- paste("chr", sapply(bins.others, function(x) strsplit(x, ";")[[1]][[1]]), sep = "")

bins.hash <- hash::hash(bins.others.coords, bins.others)

bins.new <- sapply(params.long.lst$H3K27me3$bin, function(b) AssignHash(x = b, jhash = bins.hash, null.fill = b))
bins.new2 <- sapply(pvals.long.lst$H3K27me3$bin, function(b) AssignHash(x = b, jhash = bins.hash, null.fill = b))

params.long.lst$H3K27me3$bin <- bins.new
pvals.long.lst$H3K27me3$bin <- bins.new2


# Get DE genes  -----------------------------------------------------------


jmark.ref <- "H3K27me3"
# jmark.ref <- "H3K4me3"
jmark.compare <- "H3K4me3"

params.ref <- params.long.lst[[jmark.ref]]
params.compare <- params.long.lst[[jmark.compare]]

params.merge <- left_join(params.ref, params.compare, by = c("bin", "param"))

ggplot(params.merge, aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.1) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
  theme_bw() + 
  geom_density_2d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# filter by pvals
bins.filt <- subset(pvals.long.lst[[jmark.ref]], pval < 10^-10)$bin

ggplot(params.merge %>% filter(bin %in% bins.filt), aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.1) + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
  theme_bw() + 
  geom_density_2d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Get imputed -------------------------------------------------------------


dat.imputed.lst <- lapply(jmarks, function(jmark){
  jmat <- t(log2(tm.result.lst[[jmark]]$topic %*% tm.result.lst[[jmark]]$term))
  # reorder
  cells.ordered <- dat.metas.reordered[[jmark]]$cell
  jmat <- jmat[, cells.ordered]
  if (jmark == "H3K27me3"){
    rownames(jmat) <- sapply(rownames(jmat), function(rname) AssignHash(x = rname, jhash = bins.hash, null.fill = rname))
  }
  return(jmat)
})




# Take top bins and make heatmap  -----------------------------------------


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

jbin <- "16:23221020-23231020;NM_001252505.1..St6gal1"

print(head(bins.top.filt.lst$ClusterBcells.Estimate, n = 50))
print(head(bins.top.filt.lst$ClusterEryths.Estimate, n = 50))
subset(params.ref, bin == jbin)

# check heatmap 



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


