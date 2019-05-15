# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/6-TF_analysis_downstream.R
# Plot activities in umap

library(dplyr)
library(ggplot2)
library(data.table)
library(heatmap3)

# source("~/data/scchic/m")
# source("/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis")
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

kchoose <- 75
maradir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_outputs_B6/mara_analysis_cluster_build95_B6_CorrPeakFilt.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K", kchoose, "/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-B6_H3K4me1.filt_0.99.center_TRUE_K", kchoose)
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


jmotif <- "Cebpb"

jmotif <- "Spib"

jmotif <- "Gata3"
jmotif <- "Zeb1"

jmotif <- "Tcf3"
jmotif <- "Foxo3"
jmotif <- "Meis1"

jmotif <- "Ebf1"
jmotif <- "Ikzf2"
jmotif <- "Spic"



# Zscore ------------------------------------------------------------------


jmark <- "H3K4me1"
traj.annot <- lapply(jtrajs, function(jtraj) data.frame(cell = trajs[[jmark]][[jtraj]]$cell, traj = jtraj, lambda = trajs[[jmark]][[jtraj]]$lambda)) %>%
  bind_rows()

# heatmap
# reorder columns by lambda
cell.ordering.dat <- traj.annot %>%
  filter(grepl(jmark, cell)) %>%
  group_by(traj) %>%
  arrange(traj, lambda)
cells <- cell.ordering.dat$cell

motifs.keep <- subset(mara.out$zscore, zscore >= 2)$motif
jmat <- subset(mara.out$act.mat, motif %in% motifs.keep)
colnames(jmat) <- gsub("\\.", "-", colnames(jmat))  # B6.13W1.BM.H3K4me1.3_378 -> B6-13W1-BM-H3K4me1-3_378
jmotifs <- jmat$motif
jmat$motif <- NULL

jmat <- as.matrix(t(scale(t(jmat), center=TRUE, scale=TRUE)))
# jmat <- as.matrix(scale(jmat, center=FALSE, scale=FALSE))
rownames(jmat) <- jmotifs

cells.keep <- cells[which(cells %in% colnames(jmat))]


# cut tree
K <- 7
jmeth <- "ward.D2"
# jmeth <- "centroid"
clusters <- hclust(dist(jmat[, cells.keep]), method = jmeth)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
assertthat::assert_that(length(cbPalette) == K)
clst <- cutree(clusters, k = K)
clst.dat <- data.frame(motif = names(clst), clstr = clst)

pdf("~/data/scchic/pdfs/B6_figures/TF_activities.pdf", useDingbats = FALSE)
  hm.out <- heatmap3(t(jmat[, cells.keep]), margins = c(5, 8), cexRow = 0.25, Colv = TRUE, Rowv = NA,
                     labRow = FALSE, scale = "column", revC = TRUE,
                     distfun = dist, hclustfun = hclust, method = jmeth)
  myplclust(clusters, clusters$labels, cbPalette[clst], main = jmeth)
  jmotifs <- c("Cebpb", "Spic", "Tal1", "Ebf1", "Ebf3", "Hmbox1", "Chd1", "Cdx1", "Gata3", "Bptf")
  for (jmotif in jmotifs){
    m <- PlotXYWithColor(dat.merge %>% filter(motif == jmotif), xvar = "umap1", yvar = "umap2", cname = "activity", jtitle = jmotif)
    print(m)
  }
dev.off()


# Correlate with expression  ----------------------------------------------





