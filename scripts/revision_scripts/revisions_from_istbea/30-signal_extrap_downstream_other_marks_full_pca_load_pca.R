# Jake Yeung
# Date of Creation: 2022-05-04
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/30-signal_extrap_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(velocyto.R)


# Load metadata -----------------------------------------------------------

jmarksall <- c("k4me3", "k9me3", "k27me3", "k4me1"); names(jmarksall) <- jmarksall
# jmarks <- c(jmark); names(jmarks) <- jmarks
# jmark <- "k4me3"
# jmark <- "k27me3"
# jmark <- "k9me3"
# jmark <- "k4me1"
# jmark <- args$mark
# outdir <- args$outdir

# jmark <- "k27me3"
# jmark <- "k9me3"
jmark <- "k4me3"
# jmark <- ""
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_play"
dir.create(outdir)

# outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits"

outpdf <- file.path(outdir, paste0("downstream_manual_signal_extrap_arrows.", jmark, ".", Sys.Date(), ".pdf"))
pdf(outpdf, useDingbats = FALSE)



dat.meta.lst <- lapply(jmarksall, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta)
})

# Load signal and extrap  -------------------------------------------------

inf.rds <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/signal_extrap_vectors_regions_ctypes_from_args.", jmark, ".rds")
dat.rds <- readRDS(inf.rds)

mat.signal <- lapply(dat.rds, function(jdat){
  jdat$signal.vec
}) 
mat.signal <- do.call(rbind, mat.signal)

mat.signalextrap <- lapply(dat.rds, function(jdat){
  jdat$signal.extrap.vec
})
mat.signalextrap <- do.call(rbind, mat.signalextrap)



# Do PCA ------------------------------------------------------------------

print("PCA...")
# pca.out <- prcomp(t(mat.signal), center = TRUE, scale. = TRUE)
inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/pca_output.", jmark, ".2022-04-20.rds")
pca.out <- readRDS(inf.pca.out)
# readRDS()

print("PCA... done")
saveRDS(pca.out, file = file.path(outdir, paste0("downstream_pca_out.", jmark, ".", Sys.Date(), ".rds")))

varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)

dat.pca.init <- data.frame(cell = rownames(pca.out$x), pc1 = pca.out$x[, 1], pc2 = pca.out$x[, 2], stringsAsFactors = FALSE) %>%
  left_join(., dat.meta.lst[[jmark]]) %>%
  rowwise() %>%
  mutate(status = "current")

mat.extrap.scaled <- sweep(x = mat.signalextrap, MARGIN = 1, STATS = pca.out$center, FUN = "-")
# mat.extrap.scaled <- sweep(x = mat.extrap.scaled, MARGIN = 1, STATS = pca.out$scale, FUN = "/")
mat.extrap.pca <- t(mat.extrap.scaled) %*% pca.out$rotation


dat.extrap.pca <- data.frame(cell = rownames(mat.extrap.pca), pc1 = mat.extrap.pca[, 1], pc2 = mat.extrap.pca[, 2]) %>%
  left_join(., dat.meta.lst[[jmark]]) %>%
  mutate(status = "extrap")

dat.pca.merge.wide <- left_join(dat.pca.init, subset(dat.extrap.pca, select = c(cell, pc1, pc2)), by = "cell")

ggplot(dat.pca.init, aes(x = pc1, y = pc2, color = colcode)) +
  geom_point() +
  theme_bw() +
  scale_color_identity() + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


m <- ggplot(dat.pca.merge.wide, aes(x = pc1.x, y = pc2.x, xend = pc1.y, yend = pc2.y, color = colcode)) +
  geom_point(size = 0.5, alpha = 0.35) +
  geom_segment(arrow = arrow(length=unit(0.8, "mm")), alpha = 0.35) +
  scale_color_identity() + 
  theme_bw() +
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)

subset(dat.pca.merge.wide, is.na(pc2.y))

dat.pca.merge.wide.filt <- subset(dat.pca.merge.wide, !is.na(pc1.y) & !is.na(pc2.y))

# Plot on grid  -----------------------------------------------------------

# create grid limits
# copied from (LaManno et al: show.grid.flow https://github.com/velocyto-team/velocyto.R/blob/master/R/momentum_routines.R)
grid.n  <- 25
jfactor <- 0.5
pos <- as.data.frame(subset(dat.pca.merge.wide.filt, select = c(pc1.x, pc2.x)))
ppos <- as.data.frame(subset(dat.pca.merge.wide.filt, select = c(pc1.y, pc2.y)))

# arrow estimates for each cell
ars <- data.frame(pos[,1],pos[,2],ppos[,1],ppos[,2])
colnames(ars) <- c('x0','y0','x1','y1')
arsd <- data.frame(xd=ars$x1-ars$x0,yd=ars$y1-ars$y0)
rownames(ars) <- rownames(arsd) <- rownames(pos);

rownames(pos) <- dat.pca.merge.wide.filt$cell

rx <- range(c(range(ars$x0),range(ars$x1)))
ry <- range(c(range(ars$y0),range(ars$y1)))
gx <- seq(rx[1],rx[2],length.out=grid.n)
gy <- seq(ry[1],ry[2],length.out=grid.n)

grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2

min.grid.cell.mass <- 1
garrows <- do.call(rbind,lapply(gx,function(x) {
  # cell distances (rows:cells, columns: grid points)
  cd <- sqrt(outer(pos[,2],-gy,'+')^2 + (x-pos[,1])^2)
  cw <- dnorm(cd,sd=grid.sd)
  # calculate x and y delta expectations
  gw <- Matrix::colSums(cw)
  cws <- pmax(1,Matrix::colSums(cw));
  gxd <- jfactor * Matrix::colSums(cw*arsd$xd)/cws
  gyd <- jfactor * Matrix::colSums(cw*arsd$yd)/cws
  
  al <- sqrt(gxd^2+gyd^2);
  vg <- gw>=min.grid.cell.mass & al>=min.arrow.size
  
  out <- cbind(rep(x, sum(vg)), gy[vg], x+gxd[vg], gy[vg]+gyd[vg])
}))
colnames(garrows) <- c('pc1.x','pc2.x','pc1.y','pc2.y')

dat.garrows <- data.frame(garrows)

ggplot(mapping = aes(x = 1 * pc1.x, y = pc2.x, xend = 1 * pc1.y, yend = pc2.y)) +
  geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = 0.5, size = 0.5) +
  geom_segment(arrow = arrow(length=unit(0.8, "mm")), data = dat.garrows, alpha = 0.5) +
  scale_color_identity() + 
  theme_bw() +
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

print("PCA with arrows done")
save(dat.garrows, dat.pca.merge.wide, file = file.path(outdir, paste0("downstream_pca_merged_arrows.", jmark, ".", Sys.Date(), ".RData")))


dev.off()