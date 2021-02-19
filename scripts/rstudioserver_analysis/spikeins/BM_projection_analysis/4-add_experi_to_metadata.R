# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_projection_analysis/4-add_experi_to_metadata.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new"))

assertthat::assert_that(dir.exists(indir))

outdir <- file.path(indir, "add_experi")
dir.create(outdir)

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.metadata.2020-12-27.txt")
  dat.meta.tmp <- fread(file.path(indir, fname)) %>%
    rowwise() %>%
    mutate(experi = ClipLast(x = cell, jsep = "_"))
  return(dat.meta.tmp)
})

# Write output ------------------------------------------------------------

for (jmark in jmarks){
  outfname <- paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.metadata.", Sys.Date(), ".with_experi.txt")
  outf <- file.path(outdir, outfname)
  print(outf)
  fwrite(dat.metas[[jmark]], file = outf, sep = "\t", na = "NA", quote = FALSE)
}



