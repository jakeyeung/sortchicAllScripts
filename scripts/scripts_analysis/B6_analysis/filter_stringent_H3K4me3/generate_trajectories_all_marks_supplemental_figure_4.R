# Jake Yeung
# Date of Creation: 2019-06-05
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/remake_figure_2.R
# Remake figures for main Figure 2.
# 

rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)
source("scripts/Rfunctions/PlotFunctions.R")



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- "H3K4me3"
# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.RData"
load(inf, v=T)

inf.traj.h3k4me3 <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.traj.h3k4me3, v=T)
trajs.stringent <- trajs$H3K4me3

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)
out.objs.stringent <- out.objs
tm.result.stringent <- posterior(out.objs$out.lda)

inf.lda.all <- "/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
load(inf.lda.all, v=T)


out.objs$H3K4me3 <- out.objs.stringent

inf.objs <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.objs, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent$H3K4me3
trajs$H3K4me3 <- trajs.stringent


# Plot the 4 umaps  -------------------------------------------------------

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K27me3 = "darkorange1", H3K9me3 = "red1")

ncells <- lapply(jmarks, function(jmark) return(length(unique(dat.umap.long.trajs[[jmark]]$cell))))

mlst <- lapply(jmarks, function(jmark) PlotXYNoColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", jcol = colsvec[[jmark]], jsize = 1) + ggtitle(ncells[[jmark]]))

# Plot with louvains ------------------------------------------------------

trajnames <- c("eryth", "granu", "lymphoid")
# trajname <- "granu"
# get colors
colhash <- GetTrajColors(as.hash = TRUE, add.mega = FALSE)

pdf("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/trajectories_stringent/trajectories_all_summarized.pdf", useDingbats = FALSE)
for (jmark in jmarks){
  dat.umap.sub <- dat.umap.long.trajs[[jmark]]
  m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 2) + 
    ggtitle(paste(jmark))
  for (trajname in trajnames){
    m <- m + geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = 3, color = colhash[[trajname]])
  }
  print(m)
}
dev.off()

