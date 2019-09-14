# Jake Yeung
# Date of Creation: 2019-08-15
# File: ~/projects/scchic/scripts/scripts_analysis/revisions/find_TF_regulators.R
# Find TF regulators mentioned by reviewers


library(dplyr)
library(ggplot2)
library(data.table)
library(heatmap3)
library(ggrepel)
library(topicmodels)
library(hash)
library(JFuncs)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.H3K4me1.RData"
assertthat::assert_that(file.exists(inf))
load(inf, v=T)


# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)

# Load MARA  --------------------------------------------------------------

kchoose <- 50
maradir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K", kchoose, "/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K", kchoose)
assertthat::assert_that(dir.exists(maradir))

mara.out <- LoadMARA(maradir, make.cnames = FALSE)
mara.out$act.long$cell <- gsub("\\.", "-", mara.out$act.long$cell)


dat.merge <- left_join(dat.umap.long, mara.out$act.long)

jmotif <- "Irf4"
jmotif <- "Gata1"
jmotif <- "Gata2"
jmotif <- "Klf1"

# do all Gatas
jmotifs <- unique(subset(dat.merge, grepl("Gata", motif))$motif)
jmotifs <- unique(subset(dat.merge, grepl("Klf", motif))$motif)

subset(mara.out$zscores, motif %in% jmotifs)

jmotif <- "Gata1"
jmotif <- "Klf4"

jmotif <- "Klf12"
# for (jmotif in jmotifs){
  m <- PlotXYWithColor(dat.merge %>% filter(motif == jmotif), xvar = "umap1", yvar = "umap2", cname = "activity", jtitle = jmotif, jcol = scales::muted('darkred'))
  print(m)
# }


