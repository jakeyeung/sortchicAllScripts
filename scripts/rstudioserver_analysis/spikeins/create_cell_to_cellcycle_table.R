# Jake Yeung
# Date of Creation: 2020-08-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/create_cell_to_cellcycle_table.R
# Split bam by cell cycle

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.cluster_tables"

inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")
assertthat::assert_that(file.exists(inf.spike))

load(inf.spike, v=T)

dat.annot <- AddCellCycleLabel.bydat(dat.spikeins.mat) %>%
  rowwise() %>%
  mutate(mark = strsplit(samp, split = "-")[[1]][[3]])


# Load cleaned matrix as a cell filter ------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  print(jmark)
  
  outf <- file.path(outdir, paste0(jmark, ".cell_cycle_cluster_tables.txt"))
  
  inf.mat <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, ".G1_G2_S.topn_5000.rds"))
  mat <- readRDS(inf.mat)
  
  dat.annot.sub <- subset(dat.annot, mark == jmark) %>%
    dplyr::select(c(samp, cellcycle.str)) %>%
    dplyr::rename(cluster = cellcycle.str,
                  cell = samp) %>%
    mutate(cluster = gsub(pattern = "/", replacement = "", x = cluster))
  
  fwrite(dat.annot.sub, file = outf, sep = "\t")
  
}




