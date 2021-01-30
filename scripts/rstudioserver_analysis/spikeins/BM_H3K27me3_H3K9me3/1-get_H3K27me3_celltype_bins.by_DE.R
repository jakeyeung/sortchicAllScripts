# Jake Yeung
# Date of Creation: 2021-01-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/1-get_H3K27me3_celltype_bins.by_DE.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

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


# Get celltype specific genes  --------------------------------------------

qval.cutoff <- 1e-10
params.filt <- subset(params.long, padj < qval.cutoff) %>%
  arrange(padj)
print(dim(params.filt))

bins.filt <- unique(params.filt$bin)

ggplot(params.long %>% filter(bin %in% bins.filt), aes(x = estimate, fill = param)) + 
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# print(params.long)
# jbin <- "chr1:12800000-12850000"
# jbin <- bins.filt[[2]]
jbin <- bins.filt[[1]]
print(jbin)

jsub <- subset(params.long, bin == jbin)

print(jsub)



# Get granu-specific regions ----------------------------------------------

jparams <- unique(params.long$param); names(jparams) <- jparams

# jparam <- "ClusterGranulocytes.Estimate"
params.annot.lst <- lapply(jparams, function(jparam){
  print(jparam)
  bins.filt.byctype <- subset(params.long %>% arrange(padj), param == jparam & padj < qval.cutoff)$bin
  params.long.annot <- params.long %>%
    filter(bin %in% bins.filt.byctype) %>%
    rowwise() %>%
    mutate(is.ctype = param == jparam) %>%
    group_by(bin, is.ctype) %>%
    filter(abs(log2fc) < 10) %>%
    summarise(log2fc.mean = mean(log2fc)) %>%
    group_by(bin) %>%
    summarise(log2fc.diff = log2fc.mean[[2]] - log2fc.mean[[1]]) %>%
    arrange(desc(log2fc.diff)) %>%
    mutate(param = jparam)
})


lapply(params.annot.lst, head)
lapply(params.annot.lst, tail)

# add HSPCs
jparam.toadd <- "ClusterHSPCs.Estimate"

params.sum <- params.long %>%
  group_by(bin, mark) %>%
  summarise(estimate = mean(estimate),
            se = sqrt(sum(se ^ 2)) / length(se)) %>%
  rowwise() %>%
  mutate(z = estimate / se, 
         pval.param = exp(-0.717 * z - 0.416 * z ^ 2)) %>%
  arrange(pval.param)
params.sum$padj <- p.adjust(params.sum$pval.param, method = "BH")

bins.filt.byhspcs <- subset(params.sum %>% arrange(padj), padj < qval.cutoff)$bin

params.long.annot.hspcs <- params.long %>%
  filter(bin %in% bins.filt.byhspcs) %>%
  group_by(bin) %>%
  # group_by(bin, is.ctype) %>%
  filter(abs(log2fc) < 10) %>%
  summarise(log2fc.mean = mean(log2fc)) %>%
  group_by(bin) %>%
  summarise(log2fc.diff = -1 * log2fc.mean) %>%  # switch sign positive log2fc means high in HSPCs, negative log2fc means low in HSPCs (like Hox)
  arrange(desc(log2fc.diff)) %>%
  mutate(param = jparam.toadd)


params.annot.lst[[jparam.toadd]] <- params.long.annot.hspcs

# Take top n and make heatmap  --------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"


# load metas


# load metas 

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



# Filter bins and make heatmap  -------------------------------------------

ctypes.bins <- paste("Cluster", c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs"), ".Estimate", sep = "")
names(ctypes.bins) <- ctypes.bins

assertthat::assert_that(all(names(params.annot.lst) %in% ctypes.bins))
assertthat::assert_that(all(ctypes.bins %in% names(params.annot.lst)))

# keepn <- 150
keepn <- 150
params.annot.filt.lst <- lapply(ctypes.bins, function(jparam.tmp){
  jsub.top <- params.annot.lst[[jparam.tmp]] %>%
    arrange(log2fc.diff)
    # arrange(desc(log2fc.diff))
  assertthat::assert_that(nrow(jsub.top) > keepn)
  return(jsub.top[1:keepn, ])
})


bins.byctypes <- unlist(lapply(params.annot.filt.lst, function(jsub) jsub$bin), use.names = FALSE)

jmark <- "H3K27me3"
# jmark <- "H3K9me3"
# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
colsidecolors <- dat.metas.reordered[[jmark]]$clustercol
bins.byctypes.filt <- bins.byctypes[bins.byctypes %in% rownames(dat.imputed.lst[[jmark]])]
print(length(bins.byctypes))
print(length(bins.byctypes.filt))

jmat.filt <- dat.imputed.lst[[jmark]][bins.byctypes.filt, ]

heatmap3::heatmap3(jmat.filt, labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark))


# Get raw cuts ------------------------------------------------------------

jmat.raw <- outs.all.lst[[jmark]]$count.mat[bins.byctypes.filt, ]

# normalize by the bins you chose
jmat.raw.norm.bins <- sweep(jmat.raw, MARGIN = 2, STATS = colSums(jmat.raw), FUN = "/")
jmat.raw.norm.all <- count.mat.lst[[jmark]][bins.byctypes.filt, ]

jmat.raw.norm.bins.pseudogene <- lapply(params.annot.filt.lst, function(jdat){
  unlist(colSums(jmat.raw.norm.bins[jdat$bin, ]))
  # unlist(colSums(jmat.raw.norm.all[jdat$bin, ]))
})
jmat.raw.norm.bins.pseudogene <- log2(do.call(rbind, jmat.raw.norm.bins.pseudogene) + 1)

# log and handle outliers

mat.adj <- jmat.raw.norm.bins.pseudogene
plot(density(mat.adj))

mat.adj <- t(apply(mat.adj, 1, function(jrow) DescTools::Winsorize(jrow, probs = c(0.01, 0.99))))
mat.adj <- apply(mat.adj, 2, function(jcol) DescTools::Winsorize(jcol, probs = c(0.01, 0.99)))
plot(density(mat.adj))

heatmap3::heatmap3(mat.adj, labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", revC = TRUE, main = paste0(jmark))


# Check bin ---------------------------------------------------------------


jmat.filt.long <- data.table::melt(jmat.filt)
colnames(jmat.filt.long) <- c("bin", "cell", "val")
jmat.filt.long <- left_join(jmat.filt.long, dat.metas.reordered[[jmark]])


jmat.raw.long <- data.table::melt(as.matrix(jmat.raw.norm.bins[bins.byctypes.filt, ]))
colnames(jmat.raw.long) <- c("bin", "cell", "val")
jmat.raw.long <- left_join(jmat.raw.long, dat.metas.reordered[[jmark]])


# jbins <- "chr11:44600000-44650000"
# jbins <- "chr11:44600000-44650000"
jclst <- "ClusterBcells.Estimate"
jbins <- params.annot.filt.lst[[jclst]]$bin
ggplot(jmat.filt.long %>% filter(bin %in% jbins), aes(x = ctype, y = val)) + 
  geom_boxplot() + 
  ggtitle(jclst) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmat.raw.long %>% filter(bin %in% jbins), aes(x = ctype, y = val)) + 
  scale_y_log10() +
  geom_boxplot() + 
  ggtitle(jclst) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(jmat.raw.long %>% filter(bin %in% jbins) %>% group_by(ctype) %>% summarise(frac.nonzeros = sum(ceiling(val)) / length(val)), 
       aes(x = ctype, y = frac.nonzeros)) + 
  scale_y_log10() +
  geom_boxplot() + 
  ggtitle(jclst) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check log2fc estimates


jmark2 <- "H3K9me3"
inf.fits2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark2, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
# inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.", jmark, ".2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData")
load(inf.fits2, v=T)
jfits.lst2 <- jfits.lst

params.long2 <- SummarizeParamsPvalues(jfits.lst2, jmark = jmark2, paramname = "Cluster") %>%
  mutate(log2fc = estimate / log(2))
params.long2$padj <- p.adjust(params.long2$pval.param)
jnames2 <- names(jfits.lst2); names(jnames2) <- jnames2
pvals.long2 <- lapply(jnames2, function(jname){
  x <- jfits.lst2[[jname]]
  xvec <- x$pval
  data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# merge k27me3 and k9me3 
print(unique(params.long2$param))
print(unique(params.long$param))

params.common <- intersect(unique(params.long2$param), unique(params.long$param))

params.long.merge <- left_join(params.long2, params.long %>% filter(param %in% params.common), by = c("bin", "param"))

ggplot(params.long.merge, aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.25) + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.long.merge %>% filter(padj.x < 10^-8), aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.25) + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(params.long.merge %>% filter(padj.y < 10^-30 & log2fc.y > 0), aes(x = log2fc.x, y = log2fc.y)) + 
  geom_point(alpha = 0.25) + 
  geom_density_2d() + 
  facet_wrap(~param) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()