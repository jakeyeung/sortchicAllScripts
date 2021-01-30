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



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load  -------------------------------------------------------------------

indir.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3"

params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster")
})

pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

# pvalcutoff <- 10^-10

params.wide.lst <- lapply(jmarks, function(jmark){
  params.long <- params.lst[[jmark]]
  # bins.keep <- subset(pvals.lst[[jmark]], pval < pvalcutoff)$bin
  # take top 1000 bins
  bins.keep <- (pvals.lst[[jmark]] %>% arrange(pval))$bin[1:1000]
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


# Do hierarchical clustering on each --------------------------------------


dat.metas <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned/cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  fread(inf.meta)
})

# check batch effects
mlst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~jrep)
  print(m)
})



# Filter for celtype specific genes ---------------------------------------

# 
# 
# # maybe load the robjs? 
# 
# inf.matobjs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K27me3.2021-01-08.rearranged.RData"
# load(inf.matobjs, v=T)
# 


