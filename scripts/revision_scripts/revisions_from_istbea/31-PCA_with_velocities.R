# Jake Yeung
# Date of Creation: 2022-05-05
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/31-PCA_with_velocities.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

GetArrows <- function(dat.pca.merge.wide.filt, grid.n = 40, jfactor = 1){
  # grid.n  <- 25
  # jfactor <- 0.5
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
  
  print(dim(garrows))
  dat.garrows <- data.frame(garrows)
  return(dat.garrows)
}


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/pcas_with_velocities"
outpdf <- file.path(outdir, paste0("pcas_with_velocities.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)


# Load stainings --------------------------------------------------------

# inf.stainings <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/imputation_sca1_kit_lin_no_cutoff_shift_sl4_only.2022-05-05.RData"
# load(inf.stainings, v=T)

# saveRDS(dat.sub.impute.knn.lst, file = paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat.sub.impute.knn.lst.", Sys.Date(), ".rds"))
# jcheck <- dat.sub.impute.knn.lst$k4me1

dat.sub.impute.knn.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat.sub.impute.knn.lst.2022-05-06.rds")

cell2ptime.lsk <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/cell2ptime_lsk.2022-05-07.rds")


# Load PCA with arrows  ---------------------------------------------------


jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew


jalphadots <- 1
jalphaarrows <- 0.5
jsizedots <- 1
jsizearrows <- 0.5
janglearrows <- 45
jlengtharrows <- 1.6

# print(unique(dat.pca.merge.wide$ctype.from.LL))
ctypes.ordered <- c("Granulocytes", "Bcells", "Eryths", "NKs", "Monocytes", "DCs", "Basophils", "MEP", "CMP", "pDCs", "MPPs", "LT", "ST", "HSCs")
names(ctypes.ordered) <- ctypes.ordered

# Do K27me3 ---------------------------------------------------------------

# load PCA to get varexplained

jmark <- "k27me3"


inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
pca.out <- readRDS(inf.pca.out)
inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
load(inf.objs, v=T)

dat.pca.merge.wide$ctype.from.LL <- factor(dat.pca.merge.wide$ctype.from.LL, levels = ctypes.ordered)

varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)

dat.arrows.pca <- GetArrows(dat.pca.merge.wide, grid.n = 30, jfactor = 1)
ggplot(mapping = aes(x = -1 * pc1.x, y = pc2.x, xend = -1 * pc1.y, yend = pc2.y)) +
  geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# add staining 
colors_and_legend <- tricolore::Tricolore(dat.sub.impute.knn.lst[[jmark]], 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
dat.sub.impute.knn.lst[[jmark]]$rgb <- colors_and_legend$rgb


dat.pca.merge.wide.annot <- left_join(dat.pca.merge.wide, subset(dat.sub.impute.knn.lst[[jmark]], select = c(cell, ckit_f, sca1_f, lin_f, rgb))) %>%
  rowwise() %>%
  mutate(rgb = ifelse(is.na(rgb), "grey85", rgb), 
         ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
  arrange(!is.na(sca1_f), sca1_f)


ggplot(mapping = aes(x = -1 * pc1.x, y = pc2.x, xend = -1 * pc1.y, yend = pc2.y)) +
  geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

eryth.traj <- c("HSCs", "ST", "LT", "MPPs", "MEP", "Eryths")
myeloid.traj <- c("HSCs", "ST", "LT", "MPPs", "CMP", "Basophils", "Granulocytes", "Monocytes", "DCs")
lymphoid.traj <- c("HSCs", "ST", "LT", "MPPs", "Bcells", "NKs")

sort(unique(dat.pca.merge.wide.annot$ctype.from.LL))
ggplot(subset(dat.pca.merge.wide.annot, ctype.from.LL %in% eryth.traj), aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(dat.pca.merge.wide.annot, ctype.from.LL %in% myeloid.traj), aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(dat.pca.merge.wide.annot, ctype.from.LL %in% lymphoid.traj), aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = pc2.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Do K9me3 ----------------------------------------------------------------



jmark <- "k9me3"

inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
pca.out <- readRDS(inf.pca.out)
inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
load(inf.objs, v=T)

dat.pca.merge.wide$ctype.from.LL <- factor(dat.pca.merge.wide$ctype.from.LL, levels = ctypes.ordered)

varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)

dat.arrows.pca <- GetArrows(dat.pca.merge.wide, grid.n = 30, jfactor = 1)
ggplot(mapping = aes(x = -1 * pc1.x, y = pc2.x, xend = -1 * pc1.y, yend = pc2.y)) +
  geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


colors_and_legend <- tricolore::Tricolore(dat.sub.impute.knn.lst[[jmark]], 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
dat.sub.impute.knn.lst[[jmark]]$rgb <- colors_and_legend$rgb


dat.pca.merge.wide.annot <- left_join(dat.pca.merge.wide, subset(dat.sub.impute.knn.lst[[jmark]], select = c(cell, ckit_f, sca1_f, lin_f, rgb))) %>%
  rowwise() %>%
  mutate(rgb = ifelse(is.na(rgb), "grey85", rgb), 
         ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
  arrange(!is.na(sca1_f), sca1_f)


ggplot(mapping = aes(x = -1 * pc1.x, y = pc2.x, xend = -1 * pc1.y, yend = pc2.y)) +
  geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = pc2.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





# k4me3 -------------------------------------------------------------------



jmark <- "k4me3"

inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
pca.out <- readRDS(inf.pca.out)
inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
load(inf.objs, v=T)

dat.pca.merge.wide$ctype.from.LL <- factor(dat.pca.merge.wide$ctype.from.LL, levels = ctypes.ordered)

varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)

dat.arrows.pca <- GetArrows(dat.pca.merge.wide, grid.n = 30, jfactor = 1)
ggplot(mapping = aes(x = 1 * pc1.x, y = -1 * pc2.x, xend = 1 * pc1.y, yend = -1 * pc2.y)) +
  geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

colors_and_legend <- tricolore::Tricolore(dat.sub.impute.knn.lst[[jmark]], 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
dat.sub.impute.knn.lst[[jmark]]$rgb <- colors_and_legend$rgb

dat.pca.merge.wide.annot <- left_join(dat.pca.merge.wide, subset(dat.sub.impute.knn.lst[[jmark]], select = c(cell, ckit_f, sca1_f, lin_f, rgb))) %>%
  rowwise() %>%
  mutate(rgb = ifelse(is.na(rgb), "grey85", rgb), 
         ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
  arrange(!is.na(sca1_f), sca1_f)

ggplot(mapping = aes(x = 1 * pc1.x, y = -1 * pc2.x, xend = 1 * pc1.y, yend = -1 * pc2.y)) +
  geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = 1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = -1 * pc2.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# k4me1 -------------------------------------------------------------------



jmark <- "k4me1"

inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
pca.out <- readRDS(inf.pca.out)
inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
load(inf.objs, v=T)

dat.pca.merge.wide$ctype.from.LL <- factor(dat.pca.merge.wide$ctype.from.LL, levels = ctypes.ordered)

varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)

dat.arrows.pca <- GetArrows(dat.pca.merge.wide, grid.n = 30, jfactor = 1.5)

ggplot(mapping = aes(x = -1 * pc1.x, y = 1 * pc2.x, xend = -1 * pc1.y, yend = 1 * pc2.y)) +
  geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



colors_and_legend <- tricolore::Tricolore(dat.sub.impute.knn.lst[[jmark]], 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
dat.sub.impute.knn.lst[[jmark]]$rgb <- colors_and_legend$rgb

dat.pca.merge.wide.annot <- left_join(dat.pca.merge.wide, subset(dat.sub.impute.knn.lst[[jmark]], select = c(cell, ckit_f, sca1_f, lin_f, rgb))) %>%
  rowwise() %>%
  mutate(rgb = ifelse(is.na(rgb), "grey85", rgb), 
         ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
  arrange(!is.na(sca1_f), sca1_f)

ggplot(mapping = aes(x = -1 * pc1.x, y = pc2.x, xend = -1 * pc1.y, yend = pc2.y)) +
  geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = pc2.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  ggtitle(jmark) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())







dev.off()



# Fit branch --------------------------------------------------------------

