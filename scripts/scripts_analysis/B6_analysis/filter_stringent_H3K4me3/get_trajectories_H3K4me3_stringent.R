# Jake Yeung
# Date of Creation: 2019-06-04
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/get_trajectories_H3K4me3_stringent.R
# Get trajectories from umap 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(umap)
library(igraph)
library(hash)

library(JFuncs)
library(princurve)

source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")


# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmark <- "H3K4me3"

outdir <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent"
outdirplot <- "/Users/yeung/data/scchic/pdfs/B6_figures/trajectories_stringent"
outf <- file.path(outdir, paste0("traj_objs_H3K4me3.stringent.GranuLymphFixed.Rdata"))
outfplot <- file.path(outdirplot, "traj_objs_stringent.GranuLymphFixed.pdf")
dir.create(outdir)
dir.create(outdirplot)

pdf(outfplot, useDingbats = FALSE)

# Load louvain data -------------------------------------------------------

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
assertthat::assert_that(file.exists(inf))

load(inf, v=T)
dat.umap.long$mark <- jmark


# Init traj output objects ------------------------------------------------

trajs <- vector(mode = "list", length = 4)
names(trajs) <- jmarks
trajs.objs <- vector(mode = "list", length = 4)
names(trajs.objs) <- jmarks
dat.umap.long.trajs <- vector(mode = "list", length = 4)
names(dat.umap.long.trajs) <- jmarks


# Estimate trajectories  --------------------------------------------------


# H3K4me3 Trajectories erythryoblasts  --------------------------------------------

dat.umap.sub <- dat.umap.long %>% filter(mark == jmark)

# get erythroblasts
# we want two objects, trajs and dat.umap.long.trajs 
trajname <- "eryth"

eryth.clstr <- 3
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain == eryth.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K4me3 Trajectories lymphoid  --------------------------------------------

# trajname <- "lymphoid"
trajname <- "granu"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(8, 1, 5)

dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap2 < 1.8 & (dat.umap.sub$umap2 - 2*dat.umap.sub$umap1 < 0), TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K4me3 Trajectories granulocytes  --------------------------------------------


# trajname <- "granu"
trajname <- "lymphoid"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(6, 2, 8)  # some outlier cells? 

# remove two outlier cells
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap2 > 2 & dat.umap.sub$umap1 < 3, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K4me3 Trajectories NK  --------------------------------------------

trajname <- "nk"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(7) 

dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap2 > 1, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)



# H3K4me3 Trajectories megakar?  --------------------------------------------

trajname <- "mega"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(4, 8)  # some outlier cells? 

# remove two outlier cells
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap2 < 2 & (dat.umap.sub$umap2 - 2.5 * dat.umap.sub$umap1 + 5 > 0), TRUE, FALSE)
m <- PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5) + 
  geom_abline(slope = 2.5, intercept = -5)
print(m)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# 
# # Do islands for mutual exclusive merging later ---------------------------
# 
# # megak and granu have the progenitor cluster, redo without it for merging later
# trajname <- "GranuIsland"
# 
# PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)
# 
# louv.clstr <- c(4, 7, 9)  # some outlier cells? 
# 
# # remove two outlier cells
# dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap1 > 0, TRUE, FALSE)
# PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)
# 
# # make trajectory
# trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
# trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj
# 
# # plot output 
# m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
#   geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
#   ggtitle(paste(jmark, trajname))
# print(m)
# 
# trajname <- "MegaIsland"
# PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)
# 
# louv.clstr <- c(5, 1)  # some outlier cells? 
# 
# # remove two outlier cells
# dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
# PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)
# 
# # make trajectory
# trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
# trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj
# 
# # plot output 
# m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
#   geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
#   ggtitle(paste(jmark, trajname))
# print(m)

# Save dat.umap.sub to list -----------------------------------------------

dat.umap.long.trajs[[jmark]] <- dat.umap.sub

# save to output 
dev.off()

save(dat.umap.long.trajs, trajs, trajs.objs, file = outf)

