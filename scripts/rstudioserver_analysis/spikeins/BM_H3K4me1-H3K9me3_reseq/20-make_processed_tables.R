# Jake Yeung
# Date of Creation: 2021-03-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/20-make_processed_tables.R
# Make processed tables

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jdate <- "2021-02-11"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/double_stain_outputs"
inrdata <- file.path(indir, paste0("scchix_outputs.H3K4me1-H3K9me3.", jdate, ".setseed.RData"))

load(inrdata, v=T)

hubprefix <- "/home/jyeung/hub_oudenaarden"

outf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/tables_dblstain_final/dblstain_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".txt"))

dat.umap.long.annot.pretty <- subset(dat.umap.long.annot, select = c(cell, umap1, umap2, plate, rowcoord, colcoord, batch, stypecol, louv.repress.impute, louv.act.impute, clustercol.repress, clustercol.act)) %>%
  rowwise() %>%
  dplyr::rename(sorttype = batch, 
                assigned.H3K9me3.cluster = louv.repress.impute,
                assigned.H3K4me1.cluster = louv.act.impute)
                                       
fwrite(dat.umap.long.annot.pretty, file = outf, quote = FALSE, sep = "\t", na = "NA")

