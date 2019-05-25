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

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- jmarks[["H3K4me1"]]
inmain <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.dats <- lapply(jmarks, function(jmark) file.path(inmain, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))

# Load umaps  -------------------------------------------------------------



dat.umap.long <- lapply(inf.dats, LoadUmap) %>%
  bind_rows() %>%
  dplyr::select(-repl, -techname)


# Check FACS --------------------------------------------------------------

dat.umap.sub <- dat.umap.long %>% filter(mark == jmark)

dat.facs.filt <- LoadFACSGetLoadings(jmark = jmark)
dat.merge <- left_join(dat.umap.sub, dat.facs.filt %>% dplyr::select(cell, loadings, mark))
PlotXYWithColor(dat.merge %>% filter(!is.na(loadings)), xvar = "umap1", yvar = "umap2", cname = "loadings", cont.color = TRUE, jtitle = jmark)

# Init traj output objects ------------------------------------------------

trajs <- vector(mode = "list", length = 4)
names(trajs) <- jmarks
trajs.objs <- vector(mode = "list", length = 4)
names(trajs.objs) <- jmarks


# Trajectories erythryoblasts  --------------------------------------------

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

# Trajectories granulocytes -----------------------------------------------

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



# Trajectories lymphocytes ------------------------------------------------

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


# Trajectories megakaryocyte?? --------------------------------------------


PlotXYWithColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jtitle = jmark, jsize = 5)

trajname <- "mega"
clstr.louv <- c(7)
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


# Trajectories natural killer ---------------------------------------------

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


# Do Tcells ---------------------------------------------------------------


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


# Do tcell only -----------------------------------------------------------

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



# Do neutro-island --------------------------------------------------------


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


# Save everything ---------------------------------------------------------

save(trajs, trajs.objs, dat.umap.sub, file = file.path(outdir, "traj_objs_H3K4me1.Rdata"))


# Can we project Tcells to Bcells? ----------------------------------------
# 
# trajname <- "lymphoid"
# # try to predict on cluster 8
# pred.mat <- subset(dat.umap.sub, louvain == 8, select = c(umap1, umap2, cell))
# rownames(pred.mat) <- pred.mat$cell; pred.mat$cell <- NULL
# pred.out <- project_to_curve(x = as.matrix(pred.mat), s = trajs.objs[[jmark]][[trajname]]$s)
# pred.out$lambda.norm <- pred.out$lambda / max(trajs.objs[[jmark]][[trajname]]$lambda)
# 
# dat.path <- data.frame(cell = names(pred.out$lambda), lambda.raw = pred.out$lambda, lambda = pred.out$lambda.norm, umap1 = pred.out$s[, 1], umap2 = pred.out$s[, 2])
# # add Bcells
# dat.path <- bind_rows(dat.path, trajs[[jmark]][[trajname]] %>% dplyr::select(cell, lambda.raw, lambda, umap1, umap2))
# # dat.path <- left_join(dat.path, dat.umap.sub) %>%
# #   arrange(lambda)
# 
# # add lambda to jsub
# jsub <- left_join(dat.umap.sub, dat.path) %>% mutate(lambda.raw = ifelse(is.na(lambda.raw), 0, lambda.raw))
# m.col <- PlotXYWithColor(jsub, xvar = "umap1", yvar = "umap2", jsize = 5, cname = "lambda.raw")
# print(m.col)
# 
# ggplot(as.data.frame(pred.out$s), aes(x = umap1, y = umap2)) + geom_point()
# 
# m <- PlotXYNoColor(dat.umap.sub, xvar = "umap1", yvar = "umap2", jsize = 5) + 
#   geom_path(data = dat.path %>% arrange(lambda.raw), aes(color = lambda.raw), size = 1)  + 
#   ggtitle(paste(jmark, trajname))
# print(m)

# Include neutrophil islands to granu -------------------------------------


