# Jake Yeung
# Date of Creation: 2020-12-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/3-downstream_DE_analysis_hierarchical_clustering.R
# Summarize DE anlaysis by hierarchical clusters


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# Load  -------------------------------------------------------------------

indir.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3"

inf.obj <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.H3K4me1.2020-12-12.newannot2.witherrors.TSS.RData"
load(inf.obj, v=T)

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.", jmark, ".2020-12-12.newannot2.witherrors.TSS.RData")
  } else {
    # inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.H3K4me1.2020-12-12.newannot2.witherrors.TSS.RData")
    # inf.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3_H3K27me3_rep2_rep3reseq/poisson_fit_TSS_10000.H3K27me3.2020-12-10.newannot2.rep2_rep3seq.RData"
    inf.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_TSS.H3K27me3.2020-12-12.newannot2.witherrors.TSS.RData"
  }
  # fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  # inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster")
  
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, pvals.long = pvals.long))
})

params.lst <- lapply(out.lst, function(jout) jout$params.long)
pvals.lst <- lapply(out.lst, function(jout) jout$pvals.long)


# Reassign coord to TSS  --------------------------------------------------

library(hash)
tss.orig <- unique(params.lst$H3K4me1$bin)
bins.orig <- unique(params.lst$H3K27me3$bin)

coords <- paste("chr", sapply(tss.orig, function(x) strsplit(x, ";")[[1]][[1]]), sep = "")
tss.hash <- hash::hash(coords, tss.orig)
# tss.hash <- hash::hash(tss.orig, coords)
# tss.hash.rev <- hash::invert(tss.hash)

# check bins are in tss.hash?

# jcheck <- unlist(sapply(bins.orig, function(x) AssignHash(x = x, jhash = tss.hash.rev, null.fill = NA)))
# length(which(is.na(jcheck)))  # 53 ok 

# rename
out.lst$H3K27me3$params.long$bin <- sapply(out.lst$H3K27me3$params.long$bin, function(b) AssignHash(x = b, jhash = tss.hash, null.fill = b))
out.lst$H3K27me3$pvals.long$bin <- sapply(out.lst$H3K27me3$pvals.long$bin, function(b) AssignHash(x = b, jhash = tss.hash, null.fill = b))


# pvals.renamed <- sapply(out.lst$H3K27me3$pvals.long$bin, function(b) AssignHash(x = b, jhash = tss.hash.rev, null.fill = b))
# out.lst$params.long$H3K27me3$bin <- params.renamed
# out.lst$pvals.long$H3K27me3$bin <- pvals.renamed


# Filter for celtype specific genes ---------------------------------------

# maybe load the robjs? 

inf.matobjs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K27me3.2021-01-08.rearranged.RData"
load(inf.matobjs, v=T)

jbins <- unique(rownames(mat.adj.tmp))

params.wide.lst <- lapply(jmarks, function(jmark){
  params.long <- out.lst[[jmark]]$params.long
  # bins.keep <- subset(pvals.lst[[jmark]], pval < pvalcutoff)$bin
  # take top 1000 bins
  # bins.keep <- (pvals.lst[[jmark]] %>% arrange(pval))$bin[1:1000]
  bins.keep <- jbins
  print(paste("Keeping", length(bins.keep), "bins"))
  params.long.filt <- subset(params.long, bin %in% bins.keep)
  params.wide <- as.data.frame(data.table::dcast(params.long.filt, formula = bin ~ param, value.var = "estimate"))
  params.wide$ClusterHSPCs.Estimate <- 0
  rownames(params.wide) <- params.wide$bin
  params.wide$bin <- NULL
  return(params.wide)
})


for (jmark in jmarks){
  distout <- dist(t(as.matrix(params.wide.lst[[jmark]])), method = "euclidean")
  clstout <- hclust(distout, method = "ward.D2")
  plot(clstout, main = jmark)
}


# 
# # Do hierarchical clustering on each --------------------------------------
# 
# 
# dat.metas <- lapply(jmarks, function(jmark){
#   inf.meta <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned/cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
#   fread(inf.meta)
# })
# 
# # check batch effects
# mlst <- lapply(jmarks, function(jmark){
#   m <- ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
#     geom_point() + 
#     theme_bw() + 
#     ggtitle(jmark) + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     facet_wrap(~jrep)
#   print(m)
# })



