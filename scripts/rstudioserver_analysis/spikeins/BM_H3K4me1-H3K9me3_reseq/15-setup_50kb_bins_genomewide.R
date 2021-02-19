# Jake Yeung
# Date of Creation: 2021-01-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/15-setup_dynamic_bins.from_heatmap.R
# Take k9 dynamic bins from heatmap (fewer bins)


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load dynamic bins in heatmap (low and high) ---------------------------------------------------------------

# jlow.in.k9 <- TRUE
# jkeeptop <- 150

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})

bins.keep <- names(jfits.lst.lst$H3K9me3)

dat.k9.regions <- GetBedFromCoords(bins.keep, add.chr = FALSE, strip.chr = TRUE)

# Write outputs -----------------------------------------------------------

dat.combined.regions <- dat.k9.regions
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/regions_H3K4me1_H3K9me3_dynamic_regions"
outfile <- file.path(outdir, paste0("H3K4me1_H3K9me3_50kb_genomewide.", Sys.Date(), ".txt"))
fwrite(dat.combined.regions, file = outfile, sep = "\t", col.names = FALSE)




