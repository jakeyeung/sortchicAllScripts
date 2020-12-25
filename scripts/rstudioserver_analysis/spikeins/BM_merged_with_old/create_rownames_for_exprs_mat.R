# Jake Yeung
# Date of Creation: 2020-11-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/create_rownames_for_exprs_mat.R
# Create rownames from sitecount matrix to extract reads


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load sitecount, create bedfile (no chromo) for extracting reads  --------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

N.lst <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".txt")
  dat <- fread(inf)
  return(dat)
})

beds.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dattmp <- N.lst[[jmark]]
  jchromos <- gsub("^chr", "", sapply(dattmp$Gene.ID, GetChromo))
  jstarts <- sapply(dattmp$Gene.ID, GetStart)
  jends <- sapply(dattmp$Gene.ID, GetEnd)
  bedtmp <- data.frame(Chromo = jchromos, Start = jstarts, End = jends, peakname = dattmp$Gene.ID, stringsAsFactors = FALSE)
  return(bedtmp)
})


# Outdir  -----------------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/beds.BMAllMerged2.from_peaks.from_sitecount_mat"
# dir.create(outdir)

lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("beds_from_sitecount_matrix.", jmark, ".bed")
  fwrite(beds.lst[[jmark]], col.names = FALSE, sep = "\t", quote = FALSE, file = file.path(outdir, fname))
})




