# Jake Yeung
# Date of Creation: 2019-04-09
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/SPRING_downstream.R
# Plot coordinates of spring downstream

rm(list=ls())

library(JFuncs)
library(data.table)
library(dplyr)
library(ggplot2)
library(princurve)

ReadCoords <- function(inf, jmark, scalefac.X1=1, scalefac.X2=1){
  assertthat::assert_that(file.exists(inf))
  coords <- fread(inf, col.names = c("cell_indx_0", "X1", "X2"))
  # center and scale coordinates
  coords$X1 <- scalefac.X1 * scale(coords$X1, center = TRUE, scale = TRUE)
  coords$X2 <- scalefac.X2 * scale(coords$X2, center = TRUE, scale = TRUE)
  coords$mark <- jmark
  return(coords)
}

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")

make.plots <- TRUE


# Inits -------------------------------------------------------------------

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")
trajs.spring <- list()

# Load objs ---------------------------------------------------------------

inf.objs <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
load(inf.objs, v=T)

if (make.plots){
  pdf(paste0("~/data/scchic/pdfs/build95_trajectories_spring_", Sys.Date(), ".pdf"))
}

# Load coordinates and plot  ----------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K9me3 = "red1", H3K27me3 = "darkorange1")
names(jmarks) <- jmarks
scalefacs1 <- c(1, 1, 1, -1)
scalefacs2 <- c(-1, -1, -1, -1)

infs <- lapply(jmarks, function(jmark) return(paste0("/Users/yeung/projects/SPRING_dev/datasets/coordinates/coordinates_", jmark, "_try2.txt")))
infs.cellnames <- lapply(jmarks, function(jmark) return(paste0("/Users/yeung/projects/SPRING_dev/datasets/coordinates/cell_bcs_flat_", jmark, ".txt")))


coords.dat <- mapply(function(inf, jmark, scalefac1, scalefac2) ReadCoords(inf, jmark, scalefac1, scalefac2), 
                     infs, jmarks, scalefacs1, scalefacs2, 
                     SIMPLIFY = FALSE) %>% 
  bind_rows()
coords.dat$mark <- factor(coords.dat$mark, levels = jmarks)

cellnames <- mapply(function(inf, jmark){
  cnames <- unlist(fread(inf, header = FALSE)$V1)
  # annotate with cell indx, zero based
  indx <- seq(length(cnames)) - 1
  return(data.frame(cell_indx_0 = indx, cell = cnames, mark = jmark))
}, infs.cellnames, jmarks, SIMPLIFY = FALSE) %>%
  bind_rows() 


coords.dat <- left_join(coords.dat, cellnames)

# get cellnames from cell indx

m1 <- ggplot(coords.dat, aes(x = X1, y = X2, color = mark)) + 
  geom_point(size = 0.3) +
  scale_color_manual(values = unlist(colsvec)) + 
  facet_wrap(~mark, nrow = 1) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m1)

# add louvain colors 
louvain.long <- lapply(dat.umap.long.trajs, function(jsub) jsub %>% dplyr::select(cell, louvain, mark)) %>% bind_rows()
louvain.long$mark <- factor(louvain.long$mark, levels = jmarks)
coords.dat.join <- left_join(coords.dat, louvain.long)

# plot individual marks

jmark <- "H3K4me1"
print(head(coords.dat.join %>% filter(mark == jmark)))

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(coords.dat.join %>% filter(mark == jmark), aes(x = X1, y = X2, color = louvain)) + 
    geom_point(size = 1) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(jmark)
  return(m)
})

multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)



# Check where the neutrophils are  ----------------------------------------

# highlight neutrophils

inf.facs <- lapply(jmarks, function(jmark) file.path(paste0('/Users/yeung/data/scchic/facs/', jmark, '_index_2mice_4plates.csv')))

facs.dat <- lapply(jmarks, function(jmark){
  inf <- inf.facs[[jmark]]
  dat.facs <- read.table(inf)
  # remove bad columns
  cvar <- apply(dat.facs, 2, var)
  cols.remove <- cvar == 0
  dat.facs <- dat.facs[, !cols.remove]
  
  pca.facs <- prcomp(dat.facs, center = TRUE, scale. = TRUE)
  loadings.long <- data.frame(pca.facs$rotation, facs.feature = rownames(pca.facs$rotation), stringsAsFactors = FALSE)
  samps.long <- data.frame(pca.facs$x, samp = rownames(pca.facs$x), stringsAsFactors = FALSE)
  loadings.long$mark <- jmark
  samps.long$mark <- jmark
  return(list(loadings.long = loadings.long, samps.long = samps.long))
})

samps.long <- lapply(facs.dat, function(x) return(x$samps.long)) %>%
  bind_rows()

loadings.long <- lapply(facs.dat, function(x) return(x$loadings.long)) %>%
  bind_rows()

# select neutrophils from the samp space

ggplot(samps.long, aes(x = PC1, y = PC2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~mark)

jtopn <- 20
neutros.dat <- samps.long %>%
  group_by(mark) %>%
  arrange(desc(abs(PC1))) %>%
  mutate(indx = seq(length(PC1)),
         is.neutro = ifelse(indx <= jtopn, "zNeutro", "aNotNeutro")) %>%
  rename(cell = samp) %>%
  dplyr::select(cell, mark, PC1, is.neutro)

neutros.dat.join <- left_join(coords.dat.join, neutros.dat)
neutros.dat.join <- RankOrder(neutros.dat.join, "is.neutro", out.cname = "orderrank")

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(neutros.dat.join %>% filter(mark == jmark), aes(x = X1, y = X2, color = is.neutro, shape = is.neutro, order = orderrank)) + 
    geom_point(size = 3, alpha = 0.9) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(jmark)
  return(m)
})
# print(m.lst[[1]])

multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)


# Draw trajectories on the new UMAP ? -------------------------------------

jsize <- 2
jmark <- "H3K4me1"
trajs.mark <- list()

dat.trajs <- neutros.dat.join %>%
  filter(mark == jmark) %>%
  dplyr::rename(umap1 = X1, umap2 = X2)

m.louvain <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain, shape = is.neutro, order = orderrank)) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark)
print(m.louvain)


# Granulocytes: easy ------------------------------------------------------


ctype <- "granu"
jlouv <- c(6,7,8)
dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)



# Lymphoid trajectory: average Tcell and Bcell ----------------------------



ctype <- "tcell"
jlouv <- c(10, 9, 2)
dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
  
print(m.traj)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)

ctype <- "bcell"
jlouv <- c(2, 5, 3)
dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)

# plot both
ctype <- "lymphoid"
jlouv <- c(10, 9, 2, 5, 3)
dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 
trajs.mark[[paste0(ctype, "_aggregate")]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[["bcell"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.mark[["tcell"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)

# take average of bcell and tcell trajectories

bcell.traj <- InferTrajOnUmap(dat.trajs, cname = "bcell", init.on = "umap1", return.obj = TRUE)$pc.obj
tcell.traj <- InferTrajOnUmap(dat.trajs, cname = "tcell", init.on = "umap1", return.obj = TRUE)$pc.obj

tsub <- subset(dat.trajs, tcell == TRUE, select = c(umap1, umap2, cell)) %>% arrange(desc(umap2)); rownames(tsub) <- tsub$cell; tsub$cell <- NULL; tsub <- as.matrix(tsub) 
bsub <- subset(dat.trajs, bcell == TRUE, select = c(umap1, umap2, cell)) %>% arrange(desc(umap2)); rownames(bsub) <- bsub$cell; bsub$cell <- NULL; bsub <- as.matrix(bsub)

tcells.on.btraj <- project_to_curve(tsub, bcell.traj$s)
bcells.on.ttraj <- project_to_curve(bsub, tcell.traj$s)

tcells.on.btraj <- project_to_curve(tsub, as.matrix(trajs.mark[["bcell"]] %>% dplyr::select(umap1, umap2)))
bcells.on.ttraj <- project_to_curve(bsub, as.matrix(trajs.mark[["tcell"]] %>% dplyr::select(umap1, umap2)))

tcells.on.btraj.long <- data.frame(tcells.on.btraj$s, lambda = tcells.on.btraj$lambda, cell = rownames(tcells.on.btraj$s)) %>% 
  arrange(desc(lambda)) %>% 
  mutate(lambda = lambda / max(lambda))
bcells.on.ttraj.long <- data.frame(bcells.on.ttraj$s, lambda = bcells.on.ttraj$lambda, cell = rownames(bcells.on.ttraj$s)) %>% 
  arrange(desc(lambda)) %>%
  mutate(lambda = lambda / max(lambda))

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + 
  geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = tcells.on.btraj.long, aes(x = umap1, y = umap2, color = lambda), size = 2) + 
  geom_path(data = bcells.on.ttraj.long, aes(x = umap1, y = umap2, color = lambda), size = 2)
print(m.traj)

# calculate lymphoid trajectory by taking average of the two curves, if cell covers both trajectories, then take average of the average of the two curve

# do bcells
bcell.orig <- trajs.mark[["bcell"]] %>% dplyr::select(umap1, umap2, lambda, cell) %>% dplyr::rename(umap1.orig = umap1, umap2.orig = umap2, lambda.orig = lambda)
bcell.avg <- left_join(bcell.orig, bcells.on.ttraj.long) %>%
  dplyr::rename(umap1.new = umap1,
               umap2.new = umap2,
               lambda.new = lambda) %>%
  # mutate(umap1.bcell = (umap1.new + umap1.orig) / 2,
  #        umap2.bcell = (umap2.new + umap2.orig) / 2, 
  #        lambda.bcell = (lambda.new + lambda.orig) / 2)
  mutate(umap1 = (umap1.new + umap1.orig) / 2,
         umap2 = (umap2.new + umap2.orig) / 2, 
         lambda = (lambda.new + lambda.orig) / 2,
         celltype = "bcell")

# do tcells
tcell.orig <- trajs.mark[["tcell"]] %>% dplyr::select(umap1, umap2, lambda, cell) %>% dplyr::rename(umap1.orig = umap1, umap2.orig = umap2, lambda.orig = lambda)
tcell.avg <- left_join(tcell.orig, tcells.on.btraj.long) %>%
  dplyr::rename(umap1.new = umap1,
                umap2.new = umap2,
                lambda.new = lambda) %>%
  # mutate(umap1.tcell = (umap1.new + umap1.orig) / 2,
  #        umap2.tcell = (umap2.new + umap2.orig) / 2, 
  #        lambda.tcell = (lambda.new + lambda.orig) / 2)
  mutate(umap1 = (umap1.new + umap1.orig) / 2,
         umap2 = (umap2.new + umap2.orig) / 2, 
         lambda = (lambda.new + lambda.orig) / 2,
         celltype = "tcell")

lymphoid.avg <- full_join(bcell.avg %>% 
                            dplyr::select(cell, umap1, umap2, lambda, celltype), tcell.avg %>% 
                            dplyr::select(cell, umap1, umap2, lambda, celltype)) %>%
  group_by(cell) %>% 
  summarise(umap1 = mean(umap1),
            umap2 = mean(umap2),
            lambda = mean(lambda)) %>%
  arrange(desc(umap2)) %>%
  arrange(desc(lambda)) %>%
  dplyr::select(umap1, umap2, cell, lambda)

ggplot(lymphoid.avg, aes(x = umap1, y = umap2, color = lambda)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  xlab("X1") + ylab("X2")

#  plot avg on the umap
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + 
  geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = lymphoid.avg, inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = 2) +
  # geom_path(data = lymphoid.avg %>% filter(celltype == "tcell"), inherit.aes = FALSE, aes(x = umap1, y = umap2), size = 2, color = "darkblue") +
  # geom_path(data = lymphoid.avg %>% filter(celltype == "bcell"), inherit.aes = FALSE, aes(x = umap1, y = umap2), size = 2, color = "darkred") + 
  ggtitle("Test") + 
  xlab("X1") + ylab("X2")
print(m.traj)

# add to list
trajs.mark[[ctype]] <- lymphoid.avg


# Erythyroblasts: easy ----------------------------------------------------


ctype <- "eryth"
jlouv <- c(1, 4)
dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv & dat.trajs$umap2 > 1.25
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 
trajs.mark[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)


# Plot all trajectories together ------------------------------------------

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.mark[["granu"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  geom_path(data = trajs.mark[["lymphoid"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  geom_path(data = trajs.mark[["eryth"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))  + 
  xlab("X1") + ylab("X2")
print(m.traj)



# Finished H3K4me1: add to list -------------------------------------------

trajs.spring[[jmark]] <- trajs.mark




# Do for other marks: H3K4me3 ------------------------------------------------------

jmark <- "H3K4me3"

trajs.H3K4me3 <- list()

dat.trajs <- neutros.dat.join %>%
  filter(mark == jmark) %>%
  dplyr::rename(umap1 = X1, umap2 = X2)

m.louvain <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain, shape = is.neutro, order = orderrank)) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m.louvain)


# H3K4me3 granulocytes ----------------------------------------------------

ctype <- "granu"
jlouv <- c(2, 5)

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 

trajs.H3K4me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)



# H3K4me3 on lymphoid -----------------------------------------------------


ctype <- "lymphoid"
jlouv <- c(1, 3, 6)
# filter by a line
jslope <- 3
jint <- 0.25

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv & (dat.trajs$umap2 > jslope * dat.trajs$umap1 + jint)
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  geom_abline(slope = jslope, intercept = jint) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 

trajs.H3K4me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1", flip.lambda = TRUE)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)


# H3K4me3 on erythryoids --------------------------------------------------


ctype <- "eryth"
jlouv <- 4
# filter by a line
# jslope <- 3
# jint <- 0.25

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2")
print(m)
# fit trajectory 

trajs.H3K4me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  scale_color_manual(values = cbPalette) + 
  xlab("X1") + ylab("X2")
print(m.traj)


m.all <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K4me3[["lymphoid"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.H3K4me3[["granu"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.H3K4me3[["eryth"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  ggtitle(paste(jmark, ctype)) +
  xlab("X1") + ylab("X2")
print(m.all)



# save trajectories
trajs.spring[[jmark]] <- trajs.H3K4me3


# H3K27me3 ----------------------------------------------------------------


jmark <- "H3K27me3"

trajs.H3K27me3 <- list()

dat.trajs <- neutros.dat.join %>%
  filter(mark == jmark) %>%
  dplyr::rename(umap1 = X1, umap2 = X2)

m.louvain <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain, shape = is.neutro, order = orderrank)) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2")
  ggtitle(jmark)
print(m.louvain)



# H3K27me3 myeloid --------------------------------------------------------

ctype <- "granu"
jlouv <- c(3, 4, 6, 5)

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2")
  ggtitle(jmark)
print(m)
# fit trajectory 

trajs.H3K27me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("X1") + ylab("X2") + 
  geom_path(data = trajs.H3K27me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K27me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2") + 
  scale_color_manual(values = cbPalette)
print(m.traj)


# H3K27me3 lymphoid -------------------------------------------------------

ctype <- "lymphoid"
jlouv <- c(7, 2)

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2") + 
  ggtitle(jmark)
print(m)
# fit trajectory 

trajs.H3K27me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("X1") + ylab("X2") + 
  geom_path(data = trajs.H3K27me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K27me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2") + 
  scale_color_manual(values = cbPalette)
print(m.traj)



# H3K27me3 erythryoid -----------------------------------------------------


ctype <- "eryth"
jlouv <- c(1)

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2") + 
  ggtitle(jmark)
print(m)
# fit trajectory 

trajs.H3K27me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("X1") + ylab("X2") + 
  geom_path(data = trajs.H3K27me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K27me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2") + 
  scale_color_manual(values = cbPalette)
print(m.traj)

# plot all 

m.all <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K27me3[["lymphoid"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.H3K27me3[["granu"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.H3K27me3[["eryth"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  xlab("X1") + ylab("X2") + 
  ggtitle(paste(jmark, ctype))
print(m.all)



# save trajectories 
trajs.spring[[jmark]] <- trajs.H3K27me3

# H3K9me3 -----------------------------------------------------------------

jmark <- "H3K9me3"
trajs.H3K9me3 <- list()

dat.trajs <- neutros.dat.join %>%
  filter(mark == jmark) %>%
  dplyr::rename(umap1 = X1, umap2 = X2)

m.louvain <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain, shape = is.neutro, order = orderrank)) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2") + 
  ggtitle(jmark)
print(m.louvain)


# H3K9me3 eryth --------------------------------------------------------

# lambda.filt <- 0.65
jthres <- 0.5

ctype <- "eryth"
jlouv <- c(4, 6)

hcutoff <- 1.35

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv & dat.trajs$umap2 > hcutoff
# dat.trajs[[ctype]] <- dat.trajs$umap2 > 1
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2") + 
  ggtitle(jmark)
print(m)

# fit trajectory 

jtraj <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap1")
# compress all trajectories > threshold to be threshold, then renormlize
jtraj$lambda.tmp <- sapply(jtraj$lambda, function(l) ifelse(l > jthres, jthres, l))
# renormalize
jtraj$lambda.tmp <- jtraj$lambda.tmp / max(jtraj$lambda.tmp)
jtraj$lambda <- jtraj$lambda.tmp
jtraj$lambda.tmp <- NULL

trajs.H3K9me3[[ctype]] <- jtraj


m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K9me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2")
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K9me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2") + 
  scale_color_manual(values = cbPalette)
print(m.traj)


# H3K9me3 lymphoid  -------------------------------------------------------


ctype <- "lymphoid"
jlouv <- c(2, 1, 6)
# jint <- -2.5
# jslope <- -5

vcutoff <- -0.25

dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv & (dat.trajs$umap1 < vcutoff)
# dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv & (dat.trajs$umap2 < jslope * dat.trajs$umap1 + jint)
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark) + 
  xlab("X1") + ylab("X2") + 
  geom_vline(xintercept = vcutoff)
  # geom_abline(slope = jslope, intercept = jint)
print(m)

# fit trajectory 
trajs.H3K9me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap2", flip.lambda = TRUE)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("X1") + ylab("X2") + 
  geom_path(data = trajs.H3K9me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype)) 
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K9me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2") + 
  scale_color_manual(values = cbPalette)
print(m.traj)


ctype <- "granu"
jlouv <- c(6, 5, 3)
dat.trajs[[ctype]] <- dat.trajs$louvain %in% jlouv & dat.trajs$umap1 > vcutoff
m <- ggplot(dat.trajs, aes_string(x = "umap1", y = "umap2", color = ctype, order = "orderrank")) + 
  geom_point(size = 3, alpha = 0.9) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("X1") + ylab("X2") + 
  ggtitle(jmark)
print(m)

# fit trajectory 
trajs.H3K9me3[[ctype]] <- InferTrajOnUmap(dat.trajs, cname = ctype, init.on = "umap2", flip.lambda = FALSE)
m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("X1") + ylab("X2") + 
  geom_path(data = trajs.H3K9me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + ggtitle(paste(jmark, ctype))
print(m.traj)

m.traj <- ggplot(dat.trajs, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K9me3[[ctype]], inherit.aes = FALSE, aes(x = umap1, y = umap2), size = jsize) + ggtitle(paste(jmark, ctype)) + 
  xlab("X1") + ylab("X2") + 
  scale_color_manual(values = cbPalette)
print(m.traj)

# plot all 

m.all <- ggplot(dat.trajs, aes(x = umap1, y = umap2)) + geom_point(size = 3, color = "gray50", alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_path(data = trajs.H3K9me3[["lymphoid"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.H3K9me3[["granu"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  geom_path(data = trajs.H3K9me3[["eryth"]], inherit.aes = FALSE, aes(x = umap1, y = umap2, color = lambda), size = jsize) + 
  xlab("X1") + ylab("X2") + 
  ggtitle(paste(jmark, ctype))
print(m.all)

# save trajectory
trajs.spring[[jmark]] <- trajs.H3K9me3

if (make.plots){
  dev.off()
}

# Save objects  -----------------------------------------------------------

dat.trajs.long <- neutros.dat.join

save(trajs.spring, dat.trajs.long, file = paste0("~/data/scchic/robjs/trajectory_from_spring_", Sys.Date(), ".RData"))


