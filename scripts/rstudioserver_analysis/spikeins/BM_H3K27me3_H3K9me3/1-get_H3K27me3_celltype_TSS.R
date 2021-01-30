# Jake Yeung
# Date of Creation: 2021-01-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/1-get_H3K27me3_celltype_TSS.R
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


hubprefix <- "/home/jyeung/hub_oudenaarden"

jsuffix <- "higher"  # or lower
# jsuffix <- "lower"  # or lower

make.plots <- TRUE

if (make.plots){
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K27me3_H3K9me3_analysis"
  outpdf <- file.path(outdir, paste0("heatmap_by_TSS.", Sys.Date(), ".", jsuffix, ".pdf"))
  pdf(file = outpdf, useDingbats = FALSE)
}

# jmarks <- c("H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmark1 <- "H3K27me3"
jmark2 <- "H3K9me3"

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


# Load metas --------------------------------------------------------------






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

# params.long.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   inf.fits.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.", jmark, ".2020-12-12.newannot2.witherrors.TSS.RData"))
#   # if (jmark != "H3K27me3"){
#   # } else {
#   # }
#   load(inf.fits.tmp, v=T)
#   params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
#     mutate(log2fc = estimate / log(2))
#   params.long$padj <- p.adjust(params.long$pval.param)
#   
#   
#   return(params.long)
# })
# 


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
  if (jsuffix == "lower"){
    jdat <- jdat %>%
      arrange(log2fc.diff)
  } else if (jsuffix == "higher"){
    jdat <- jdat %>%
      arrange(desc(log2fc.diff))
  } else {
    print(paste("suffix msut be lower or higher", jsuffix))
  }
  jdat$bin[1:keepn]
})

# add HSPCs? 
params.hspcs <- params.ref %>%
  group_by(bin) %>%
  filter(abs(estimate) < 5) %>%
  summarise(log2fc.diff = mean(estimate), 
            se.diff = sqrt(sum(se ^ 2)))

if (jsuffix == "lower"){
  params.hspcs.ordered <- params.hspcs %>%
    arrange(desc(log2fc.diff))
} else if (jsuffix == "higher"){
  params.hspcs.ordered <- params.hspcs %>%
    arrange(log2fc.diff)
} else {
  print(paste("suffix msut be lower or higher", jsuffix))
}
bins.hspcs <- params.hspcs.ordered$bin[1:keepn]

# add to list
jparam.hspcs <- "ClusterHSPCs.Estimate"
bins.top.filt.lst[[jparam.hspcs]] <- bins.hspcs

jbin <- "16:23221020-23231020;NM_001252505.1..St6gal1"

print(head(bins.top.filt.lst$ClusterBcells.Estimate, n = 50))
  print(head(bins.top.filt.lst$ClusterEryths.Estimate, n = 50))
subset(params.ref, bin == jbin)

# check heatmap 


# Plot heatmap ------------------------------------------------------------


clst2col <- hash::hash(dat.metas[[jmark.ref]]$cluster, dat.metas[[jmark.ref]]$clustercol)

# ctypes.bins <- paste("Cluster", c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs"), ".Estimate", sep = "")
# ctypes.bins <- paste("Cluster", c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils"), ".Estimate", sep = "")
ctypes.raw <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs"); names(ctypes.raw) <- ctypes.raw
ctypes.bins <- paste("Cluster", ctypes.raw, ".Estimate", sep = ""); names(ctypes.bins) <- ctypes.bins
ctypes.cols <- sapply(ctypes.raw, function(x) AssignHash(x = x, jhash = clst2col, null.fill = x))
names(ctypes.cols) <- names(ctypes.bins)

# names(ctypes.bins) <- ctypes.bins
bins.keep <- unlist(bins.top.filt.lst[ctypes.bins])


colsidecolors <- dat.metas.reordered[[jmark.ref]]$clustercol

# add rowsidecols
rowsidecolors <- lapply(ctypes.bins, function(x){
  jbins.tmp <- bins.top.filt.lst[[x]]  # vector many bins
  jcols.tmp <- ctypes.cols[[x]]  # one color
  # print(jcols.tmp)
  jcols.tmp.vec <- rep(jcols.tmp, length(jbins.tmp))
}) %>%
  unlist()

heatmap3::heatmap3(dat.imputed.lst[[jmark.ref]][bins.keep, ], labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, RowSideColors = rowsidecolors, scale = "row", revC = TRUE, main = paste0(jmark.ref), ColSideLabs = "cell", RowSideLabs = "TSS")


# show another mark

jmark.tmp <- "H3K4me1"
jmark.tmp <- "H3K4me1"
jmark.tmp <- "H3K4me3"
jmark.tmp <- "H3K9me3"

for (jmark.tmp in jmarks[!jmarks %in% jmark.ref]){
  print(jmark.tmp)
  bins.keep.remove <- !bins.keep %in% rownames(dat.imputed.lst[[jmark.tmp]])
  bins.keep2 <- bins.keep[!bins.keep.remove]
  rowsidecolors2 <- rowsidecolors[!bins.keep.remove]
  heatmap3::heatmap3(dat.imputed.lst[[jmark.tmp]][bins.keep2, ], labCol = "", Rowv = NA, Colv = NA, ColSideColors = dat.metas.reordered[[jmark.tmp]]$clustercol, RowSideColors = rowsidecolors2, scale = "row", revC = TRUE, main = paste0(jmark.tmp), ColSideLabs = "cell", RowSideLabs = "TSS")
}


# 
# # Compare DEs -------------------------------------------------------------
# 
# params.ref <- params.long.lst[[jmark.ref]]
# params.compare <- params.long.lst[[jmark.tmp]]
# 
# params.merge <- left_join(params.ref, params.compare, by = c("bin", "param"))
# 
# ggplot(params.merge %>% filter(bin %in% bins.keep), aes(x = log2fc.x, y = log2fc.y)) + 
#   geom_point(alpha = 0.1) +
#   ggtitle(paste(jmark.ref, jmark.tmp)) + 
#   facet_wrap(~param) + 
#   coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
#   geom_density_2d() + 
#   theme_bw() +
#   geom_hline(yintercept = 0, linetype = "dotted") + 
#   geom_vline(xintercept = 0, linetype = "dotted") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


if (make.plots){
  dev.off()
}

