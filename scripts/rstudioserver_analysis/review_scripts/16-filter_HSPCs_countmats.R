# Jake Yeung
# Date of Creation: 2021-09-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/16-filter_HSPCs_countmats.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks.nok9 <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks.nok9) <- jmarks.nok9

# Load metadata -----------------------------------------------------------

infs.meta <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells/metadata_batch_corrected.arranged_by_lineage.shuffled.", jmark, ".2021-02-19.txt"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

dats.meta <- lapply(infs.meta, fread)

# filter HSPCs
dats.meta <- lapply(dats.meta, function(jdat){
  subset(jdat, cluster == "HSPCs")
})

# Load TSS  ---------------------------------------------------------------

# no K9me3

infs.tss.fits <- lapply(jmarks.nok9, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_TSS_10000.", jmark, ".2020-11-14.newannot2.RData"))
  if (jmark == "H3K27me3"){
    inf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3_H3K27me3_rep2_rep3reseq/poisson_fit_TSS_10000.H3K27me3.2021-02-02.newannot2.rep2_rep3seq.with_se.RData")
  }
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

mats.tss <- lapply(infs.tss.fits, function(jinf){
  load(jinf, v=T)
  return(jmat.mark)
})

mats.filt.tss <- lapply(jmarks.nok9, function(jmark){
  cells.keep <- dats.meta[[jmark]]$cell
  print(length(cells.keep))
  # cols.keep <- colnames(mats.tss[[jmark]]) %in% cells.keep
  # mats.filt.tmp <- mats.tss[[jmark]][, cols.keep]
  mats.filt.tmp <- mats.tss[[jmark]][, cells.keep]
  print(dim(mats.filt.tmp))
  return(mats.filt.tmp)
})

# Load bins  --------------------------------------------------------------

# with K9me3
infs.bins.fits <- lapply(jmarks, function(jmark){
  # if (jmark == "H3K27me3"){
  #   inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.HSPCs_vs_nonHSPCs/poisson_fit_bins.", jmark, ".2021-02-04.newannot2.witherrors.MoreBins.HSPCs_vs_nonHSPCs.total.RData"))
  # } else {
  #   inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_bins.", jmark, ".2020-12-12.newannot2.witherrors.MoreBins.RData"))
  # }
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.", jmark, ".2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

mats.bins <- lapply(infs.bins.fits, function(jinf){
  load(jinf, v=T)
  return(jmat.mark)
})

lapply(mats.bins, dim)

mats.filt.bins <- lapply(jmarks, function(jmark){
  cells.keep <- dats.meta[[jmark]]$cell
  print(length(cells.keep))
  mats.filt.tmp <- mats.bins[[jmark]][, cells.keep]
  print(dim(mats.filt.tmp))
  return(mats.filt.tmp)
})


# Write outputs -----------------------------------------------------------

# write TSS
outdir.tss <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/post_submission/HSPCs_filt_count_mats/TSS_10kb_no_k9")
outdir.bins <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/post_submission/HSPCs_filt_count_mats/genomewide_bins_50kb")

for (jmark in jmarks){
  print(jmark)
  outf.tss <- file.path(outdir.tss, paste0("HSPCs_filt_countmat_TSS_10kb.", jmark, ".", Sys.Date(), ".rds"))
  outf.bins <- file.path(outdir.bins, paste0("HSPCs_filt_countmat_bins_50kb.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(mats.filt.tss[[jmark]], outf.tss)
  saveRDS(mats.filt.bins[[jmark]], outf.bins)
}


