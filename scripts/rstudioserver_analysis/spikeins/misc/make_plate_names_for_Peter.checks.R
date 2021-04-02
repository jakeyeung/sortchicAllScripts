# Jake Yeung
# Date of Creation: 2021-02-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/misc/make_plate_names_for_Peter.R
# Make plate names for Peter

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)


hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

dat.metas <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt"))
  fread(inf)
})

dat.metas.all <-lapply(jmarks, function(jmark){
  dat.metas[[jmark]] %>%
    rowwise() %>%
    mutate(experi = ClipLast(x = cell, jsep = "_"))
})

dat.metas.filt <- lapply(dat.metas, function(jdat){
  subset(jdat, jrep != "rep1old") %>%
    rowwise() %>%
    mutate(experi = ClipLast(x = cell, jsep = "_"))
})


dat.metas.filt.old <- lapply(jmarks, function(jmark){
  if (jmark == "H3K27me3"){
    return(NULL)
  }
  jdat <- dat.metas[[jmark]]
  subset(jdat, jrep == "rep1old") %>%
    rowwise() %>%
    mutate(experi = ClipLast(x = cell, jsep = "_"))
})



dat.platenames.rep2rep3 <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/plate_names/plate_names_rep2rep3.", jmark, ".txt")
  fread(inf)
})
