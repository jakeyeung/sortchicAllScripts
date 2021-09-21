# Jake Yeung
# Date of Creation: 2021-04-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/1-write_dynamic_bins_by_mark_and_gene_regions_DE_final_table.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.dynamic_bins_TSS_TES_regions")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
bsize <- 50000

# Load genomic regions  ---------------------------------------------------

de.bins.dat.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  dat.de.bins.tmp <- fread(inf.tmp)
  jbins <- dat.de.bins.tmp$CoordOriginal
  return(jbins)
})

# make beds
de.bins.bed.lst <- lapply(jmarks, function(jmark){
  jcoords <- de.bins.dat.lst[[jmark]]
  GetBedFromCoords(jcoords, add.chr = FALSE, strip.chr = TRUE)
})


for (jmark in jmarks){
  outftmp <- file.path(outdir, paste0("dynamic_bins.50kb.corrected_DE_tables.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(de.bins.bed.lst[[jmark]], file = outftmp, sep = "\t", col.names = FALSE)
}


