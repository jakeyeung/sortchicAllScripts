# Jake Yeung
# Date of Creation: 2019-04-02
# File: ~/data/scchic/build95_primetime/pseudotime_analysis.R
# Load files for pseudotime

rm(list=ls())

library(dplyr)
library(ggplot2)
library(princurve)

source("scripts/Rfunctions/TrajFunctions.R")

# Load file  --------------------------------------------------------------


# inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.RData"
inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"

load(inf, v=T)

pdf("~/data/scchic/pdfs/build95_trajectories.pdf", useDingbats = FALSE)

# Constants ---------------------------------------------------------------

jsize <- 2

# Init objects ------------------------------------------------------------

dat.umap.long.trajs <- list()
trajs <- list()  # one for each mark

# Do pseudotime on H3K4me1 ------------------------------------------------

jmark <- "H3K4me1"

# init for H3K4me1
dat.umap.long <- dat.umap.long.new.lst[[jmark]] %>% mutate(umap2 = -umap2)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(jmark)


trajs.mark <- list()


# Do principal curves: H3K4me1 ------------------------------------------------------




# Erythropoiesis ----------------------------------------------------------

# do erythryopoesis 

ctype <- "eryth"

cells.keep <- subset(dat.umap.long, louvain %in% c(1, 4) & umap2 > 2.5)$cell
# remove outlier cells
# cells.remove <- subset(dat.umap.long, louvain %in% c(1, 4))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2")


# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)  + ggtitle(paste(jmark, ctype))

# Neutrophils -------------------------------------------------------------

ctype <- "granu"
cells.keep <- subset(dat.umap.long, (louvain %in% c(7, 8)) | (louvain == 6 & umap2 > -0.15 & umap2 < 2.5))$cell
# cells.keep <- subset(dat.umap.long, (louvain %in% c(7, 8)))$cell
# cells.keep <- subset(dat.umap.long, (louvain == 6 & umap2 > 0))$cell
# remove outlier cells
# cells.remove <- subset(dat.umap.long, louvain %in% c(1, 4))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2")

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))


# Megakaryo? --------------------------------------------------------------


ctype <- "mega"
cells.keep <- subset(dat.umap.long, (louvain == 6 & umap2 < 2.5 & umap1 > -0.1))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = TRUE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))
# Bcell -------------------------------------------------------------------

ctype <- "bcell"
cells.keep <- subset(dat.umap.long, (louvain %in% c(2, 5, 3) & umap1 < 0))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2")

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))


# Tcell -------------------------------------------------------------------

ctype <- "tcell"

cells.keep <- subset(dat.umap.long, (louvain %in% c(2, 9) & umap1 < 0))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = TRUE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))

# NK cells ----------------------------------------------------------------


ctype <- "nkcell"

cells.keep <- subset(dat.umap.long, (louvain %in% c(10) & umap1 < 0))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2")

# plot it
jsize <- 2
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))


# Summarize all trajs for H3K4me1 -----------------------------------------

m.all <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, "H3K4me1 trajectories"))

for (ctype in names(trajs.mark)){
  m.all <- m.all + geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)
}

print(m.all)

# plot but dont save yet, do lymphoid anad myeloid below

# Do lymphoid and myeloid in addition  ------------------------------------

# myeloid

ctype <- "myeloid"

cells.keep <- subset(dat.umap.long, (louvain %in% c(6, 7, 8)) & umap2 < 2.5 & umap1 > 0)$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap1", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))

# lymphoid
ctype <- "lymphoid"

cells.keep <- subset(dat.umap.long, (louvain %in% c(2, 5, 9, 3) & umap1 < 0))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap1", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))



# Plot the myeloid and lymphoid trajectories ------------------------------

m.all <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, "Eryth, Myeloid, Lymphoid trajectories"))

for (ctype in c("eryth", "myeloid", "lymphoid")){
  m.all <- m.all + geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)
}

print(m.all)

dat.umap.long.trajs[[jmark]] <- dat.umap.long

trajs[[jmark]] <- trajs.mark

# now we save


dat.umap.long.trajs[[jmark]] <- dat.umap.long

trajs[[jmark]] <- trajs.mark



# Do principal curves: H3K4me3 ------------------------------------------------------

jmark <- "H3K4me3"
dat.umap.long <- dat.umap.long.new.lst[[jmark]] %>% mutate(umap2 = -umap2)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark)


trajs.mark <- list()

# Do eryth ----------------------------------------------------------------


ctype <- "eryth"

cells.keep <- subset(dat.umap.long, louvain == 4 & umap2 > 4)$cell
# remove outlier cells
# cells.remove <- subset(dat.umap.long, louvain %in% c(1, 4))$cell
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))


# Do lymphoid lineage ----------------------------------------------------------------


ctype <- "lymphoid"
cells.keep <- subset(dat.umap.long, louvain %in% c(6, 3, 1) & umap1 < 0.5 & umap2 < 4)$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = TRUE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)  + ggtitle(paste(jmark, ctype))


# Do myeloid lineage ------------------------------------------------------


ctype <- "myeloid"

cells.keep <- subset(dat.umap.long, louvain %in% c(2, 5) & umap1 > -1)$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = TRUE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)  + ggtitle(paste(jmark, ctype))


# Summarize all trajs for H3K4me3 -----------------------------------------


m.all <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, "All Trajs"))

for (ctype in names(trajs.mark)){
  m.all <- m.all + geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)
}

print(m.all)

# rewrite the dat.umap to include celltypes
# dat.umap.long.new.lst[[jmark]] <- dat.umap.long
dat.umap.long.trajs[[jmark]] <- dat.umap.long

trajs[[jmark]] <- trajs.mark


# Do H3K27me3 -------------------------------------------------------------

jmark <- "H3K27me3"
dat.umap.long <- dat.umap.long.new.lst[[jmark]]
trajs.mark <- list()

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark))


# Myeloid ----------------------------------------------------------------


ctype <- "myeloid"

cells.keep <- subset(dat.umap.long, louvain %in% c(3, 4, 6, 5))$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)  + ggtitle(paste(jmark, ctype))


# Lymphoid ----------------------------------------------------------------


ctype <- "lymphoid"

cells.keep <- subset(dat.umap.long, louvain %in% c(7, 2))$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)  + ggtitle(paste(jmark, ctype))


# Eryth -------------------------------------------------------------------


ctype <- "eryth"

cells.keep <- subset(dat.umap.long, louvain %in% c(1) & umap2 > 5)$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)  + ggtitle(paste(jmark, ctype))



# Summarize all trajs for H3K27me3 -----------------------------------------


m.all <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + ggtitle(paste(jmark, "All Trajs"))

for (ctype in names(trajs.mark)){
  m.all <- m.all + geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)
}

print(m.all)

dat.umap.long.trajs[[jmark]] <- dat.umap.long

trajs[[jmark]] <- trajs.mark


# Init for H3K9em3 --------------------------------------------------------


jmark <- "H3K9me3"
dat.umap.long <- dat.umap.long.new.lst[[jmark]]
trajs.mark <- list()

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))


# Myeloid ------------------------------------------------------------------



ctype <- "myeloid"

cells.keep <- subset(dat.umap.long, louvain %in% c(6, 5, 3) & umap1 < 3)$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = TRUE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))



# Lymphoid ----------------------------------------------------------------


ctype <- "lymphoid"

cells.keep <- subset(dat.umap.long, louvain %in% c(6, 2, 1))$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = TRUE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))



# Eryth -------------------------------------------------------------------

ctype <- "eryth"

cells.keep <- subset(dat.umap.long, louvain %in% c(4) & umap1 > 3)$cell
# remove outlier cells
dat.umap.long[[ctype]] <- sapply(dat.umap.long$cell, function(x) x %in% cells.keep)

# label these on 

# check
m1 <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = ctype)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, ctype))
print(m1)

# filter out cells and run principal curve
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.umap.long, cname = ctype, init.on = "umap2", flip.lambda = FALSE)

# plot it
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))


# Summarize all H3K9me3  ----------------------------------------------------------


m.all <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ ggtitle(paste(jmark, "All Trajs"))

for (ctype in names(trajs.mark)){
  m.all <- m.all + geom_path(data = trajs.mark[[ctype]], aes(color = lambda), size = jsize)
}

print(m.all)

dat.umap.long.trajs[[jmark]] <- dat.umap.long

trajs[[jmark]] <- trajs.mark


# Save everything for AvO -------------------------------------------------

dev.off()

save(dat.umap.long.trajs, trajs, mara.outs, tm.result.lst, count.mat.lst, annots.lst, custom.settings.new.lst, file = paste0("/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.", Sys.Date(), ".RData"))
