# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/7-get_trajectories.R
# Get trajectories for the 4 marks 

library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(JFuncs)
library(princurve)

source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")

# Constants ---------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/robjs/B6_objs"
outdirplot <- "/Users/yeung/data/scchic/pdfs/B6_figures/trajectories"
outf <- file.path(outdir, "traj_objs_all_marks.Rdata")
outfplot <- file.path(outdirplot, "traj_objs_all_marks.pdf")

pdf(outfplot, useDingbats = FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- jmarks[["H3K4me1"]]
inmain <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.dats <- lapply(jmarks, function(jmark) file.path(inmain, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))

# Load umaps  -------------------------------------------------------------

dat.umap.long <- lapply(inf.dats, LoadUmap) %>%
  bind_rows() %>%
  dplyr::select(-repl, -techname)

# Init traj output objects ------------------------------------------------

trajs <- vector(mode = "list", length = 4)
names(trajs) <- jmarks
trajs.objs <- vector(mode = "list", length = 4)
names(trajs.objs) <- jmarks
dat.umap.long.trajs <- vector(mode = "list", length = 4)
names(dat.umap.long.trajs) <- jmarks


# H3K4me1 Check FACS --------------------------------------------------------------

dat.umap.sub <- dat.umap.long %>% filter(mark == jmark)
dat.facs.filt <- LoadFACSGetLoadings(jmark = jmark)
dat.merge <- left_join(dat.umap.sub, dat.facs.filt %>% dplyr::select(cell, loadings, mark))
PlotXYWithColor(dat.merge %>% filter(!is.na(loadings)), xvar = "umap1", yvar = "umap2", cname = "loadings", cont.color = TRUE, jtitle = jmark)

# H3K4me1 Trajectories erythryoblasts  --------------------------------------------

# H3K4me1

# check FACS
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# get erythroblasts
# we want two objects, trajs and dat.umap.long.trajs 
trajname <- "eryth"

eryth.clstr <- 2
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain == eryth.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K4me1 Trajectories granulocytes -----------------------------------------------

# H3K4me1
# check FACS
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)


trajname <- "granu"

# granu.louv <- c(1, 6, 10, 4, 5)
granu.louv <- c(1, 10, 4, 5)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% granu.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)



# H3K4me1 Trajectories lymphocytes ------------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "lymphoid"

clstr.louv <- c(3, 11)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# # try to predict on cluster 8
# pred.mat <- subset(dat.umap.sub, louvain == 8, select = c(umap1, umap2, cell))
# rownames(pred.mat) <- pred.mat$cell; pred.mat$cell <- NULL
# pred.out <- project_to_curve(x = as.matrix(pred.mat), s = trajs.objs[[jmark]][[trajname]]$s)
# pred.out$lambda.norm <- pred.out$lambda / max(trajs.objs[[jmark]][[trajname]]$lambda)


# H3K4me1 Trajectories megakaryocyte?? --------------------------------------------


PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "mega"
clstr.louv <- c(7, 10)  
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv & dat.umap.sub$umap2 > -1, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

print("DOING MEGA WITH CENTER LOCATION")
# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# H3K4me1 do mega as island -----------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "MegaIsland"
clstr.louv <- c(7)  
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

print("DOING MEGA AS ISLAND")
# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)



# H3K4me1 Trajectories natural killer ---------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "nk"
clstr.louv <- c(9)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# H3K4me1 Do Tcells ---------------------------------------------------------------


PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "Tcell"
clstr.louv <- c(3, 8)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# H3K4me1 Do tcell only -----------------------------------------------------------

# helps with clustering bams

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "TcellIsland"
clstr.louv <- c(8)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], aes(color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)



# H3K4me1 Do neutro-island --------------------------------------------------------


PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "NeutroIsland"
# clstr.louv <- c(10, 6)
clstr.louv <- c(6)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% clstr.louv, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE, flip.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE, flip.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 1)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# Save dat.umap.sub to list -----------------------------------------------

dat.umap.long.trajs[[jmark]] <- dat.umap.sub



# Initialize H3K4me3 object -----------------------------------------------

jmark <- "H3K4me3"
dat.umap.sub <- dat.umap.long %>% filter(mark == jmark)
dat.facs.filt <- LoadFACSGetLoadings(jmark = jmark)
dat.merge <- left_join(dat.umap.sub, dat.facs.filt %>% dplyr::select(cell, loadings, mark))
PlotXYWithColor(dat.merge %>% filter(!is.na(loadings)), xvar = "umap1", yvar = "umap2", cname = "loadings", cont.color = TRUE, jtitle = jmark, jsize = 5)


# H3K4me3 Trajectories erythryoblasts  --------------------------------------------

# get erythroblasts
# we want two objects, trajs and dat.umap.long.trajs 
trajname <- "eryth"

eryth.clstr <- 2
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

trajname <- "lymphoid"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(3, 8, 1)

dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
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


trajname <- "granu"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(4, 7, 9, 1)  # some outlier cells? 

# remove two outlier cells
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap1 > 0, TRUE, FALSE)
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

louv.clstr <- c(6) 

dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
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

louv.clstr <- c(5, 1)  # some outlier cells? 

# remove two outlier cells
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# Do islands for mutual exclusive merging later ---------------------------

# megak and granu have the progenitor cluster, redo without it for merging later
trajname <- "GranuIsland"

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(4, 7, 9)  # some outlier cells? 

# remove two outlier cells
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap1 > 0, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

trajname <- "MegaIsland"
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

louv.clstr <- c(5, 1)  # some outlier cells? 

# remove two outlier cells
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap2", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# Save dat.umap.sub to list -----------------------------------------------

dat.umap.long.trajs[[jmark]] <- dat.umap.sub

# Initialize H3K27me3 object -----------------------------------------------

jmark <- "H3K27me3"
dat.umap.sub <- dat.umap.long %>% filter(mark == jmark)
dat.facs.filt <- LoadFACSGetLoadings(jmark = jmark)
dat.merge <- left_join(dat.umap.sub, dat.facs.filt %>% dplyr::select(cell, loadings, mark))
PlotXYWithColor(dat.merge %>% filter(!is.na(loadings)), xvar = "umap1", yvar = "umap2", cname = "loadings", cont.color = TRUE, jtitle = jmark, jsize = 5)


# H3K27me3 Trajectories erythryoblasts  --------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)


# get erythroblasts
# we want two objects, trajs and dat.umap.long.trajs 
trajname <- "eryth"

eryth.clstr <- 3
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain == eryth.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K27me3 Trajectories lymphoid  --------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "bcell"

louv.clstr <- c(5, 7)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# H3K27me3 granu -----------------------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "granu"

louv.clstr <- c(2, 8)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)


# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# H3K27me3 mega? -----------------------------------------------------------


PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "mega"  # remove some outliers

louv.clstr <- c(1, 6)

# remove 2 outliers
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & abs(dat.umap.sub$umap2) < 4, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# Tcell traj? --------------------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "tcell"  # remove some outliers

louv.clstr <- c(4)

dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# Tcell Bcell merge -------------------------------------------------------


PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# trajname <- "TcellBcellMerge"  # remove some outliers
trajname <- "lymphoid"  # remove some outliers

louv.clstr <- c(4, 7, 5)

dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# Write to list -----------------------------------------------------------

dat.umap.long.trajs[[jmark]] <- dat.umap.sub

# Initialize H3K9me3 object -----------------------------------------------

jmark <- "H3K9me3"
dat.umap.sub <- dat.umap.long %>% filter(mark == jmark)
dat.facs.filt <- LoadFACSGetLoadings(jmark = jmark)
dat.merge <- left_join(dat.umap.sub, dat.facs.filt %>% dplyr::select(cell, loadings, mark))
PlotXYWithColor(dat.merge %>% filter(!is.na(loadings)), xvar = "umap1", yvar = "umap2", cname = "loadings", cont.color = TRUE, jtitle = jmark, jsize = 5)


# H3K9me3 Trajectories erythryblasts  --------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "eryth"

louv.clstr <- c(7)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K9me3 Trajectories lymphoid  --------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "lymphoid"

louv.clstr <- c(5, 3)
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

# H3K9me3 Trajectories granu  --------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "granu"

louv.clstr <- c(6, 1, 4)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & dat.umap.sub$umap2 < 1, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE, flip.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE, flip.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)


# H3K9me3 Trajectories megakar?  --------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "mega"

louv.clstr <- c(2)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE, flip.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE, flip.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

# H3K9me3 merge granu -----------------------------------------------------

PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "GranuMegaMerge"

louv.clstr <- c(2, 6, 1, 4)
dat.umap.sub[[trajname]] <- ifelse(dat.umap.sub$louvain %in% louv.clstr & louv.clstr & dat.umap.sub$umap2 < 1, TRUE, FALSE)
PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = trajname, cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

# make trajectory
trajs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = FALSE, get.raw.lambda = TRUE, flip.lambda = TRUE)
trajs.objs[[jmark]][[trajname]] <- InferTrajOnUmap(dat = dat.umap.sub, cname = trajname, init.on = "umap1", return.obj = TRUE, get.raw.lambda = TRUE, flip.lambda = TRUE)$pc.obj

# plot output 
m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
  geom_path(data = trajs[[jmark]][[trajname]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2)  + 
  ggtitle(paste(jmark, trajname))
print(m)

dev.off()

# Write to list -----------------------------------------------------------

dat.umap.long.trajs[[jmark]] <- dat.umap.sub


# Check and write to file  ------------------------------------------------

save(dat.umap.long.trajs, trajs, trajs.objs, file = outf)
