# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/6-TF_analysis_downstream.R
# Plot activities in umap

library(dplyr)
library(ggplot2)
library(data.table)

# source("~/data/scchic/m")
# source("/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis")
source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.H3K4me1.RData"

assertthat::assert_that(file.exists(inf))

load(inf, v=T)


# Load MARA  --------------------------------------------------------------

maradir <- "/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K50"
assertthat::assert_that(dir.exists(maradir))

mara.out <- LoadMARA(maradir, make.cnames = FALSE)
mara.out$act.long$cell <- gsub("\\.", "-", mara.out$act.long$cell)


# Merge and plot ----------------------------------------------------------

dat.merge <- left_join(dat.umap.long, mara.out$act.long)

jmotif <- "Irf4"
jmotif <- "Zeb1"
jmotif <- "Tcf3"
jmotif <- "Stat2"
jmotif <- "Yy1"
jmotif <- "Hoxa5"
jmotif <- "Rara"
jmotif <- "Sox6"
jmotif <- "Zeb1"


jmotif <- "Tal1"
jmotif <- "Cebpb"

jmotif <- "Spib"
jmotif <- "Ebf1"

PlotXYWithColor(dat.merge %>% filter(motif == jmotif), xvar = "umap1", yvar = "umap2", cname = "activity", jtitle = jmotif)


