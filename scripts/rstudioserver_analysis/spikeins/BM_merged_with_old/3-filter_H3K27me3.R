# Jake Yeung
# Date of Creation: 2020-11-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/3-filter_H3K27me3.R
# rerun LDA

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/count_mat_old_merged_with_new.H3K27me3.rds")
outf <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt/count_mat_old_merged_with_new.H3K27me3.rds")
assertthat::assert_that(file.exists(inf))

count.mat <- readRDS(inf)


# Load annots -------------------------------------------------------------

inf.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2/spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.2020-11-01.WithRelLevels.mark_H3K27me3.cell_cluster_tables.txt"
dat.annot <- fread(inf.annot)



# Remove bad cells --------------------------------------------------------

dat.annot.filt <- subset(dat.annot, !grepl("louvain", cluster))

cells.keep <- dat.annot.filt$cell

count.mat.filt <- count.mat[, cells.keep]

saveRDS(object = count.mat.filt, file = outf)

