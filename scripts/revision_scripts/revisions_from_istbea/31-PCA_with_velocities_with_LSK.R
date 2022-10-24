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
outpdf <- file.path(outdir, paste0("pcas_with_velocities_with_LSK.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)


# Load stainings --------------------------------------------------------

# inf.stainings <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/imputation_sca1_kit_lin_no_cutoff_shift_sl4_only.2022-05-05.RData"
# load(inf.stainings, v=T)

# saveRDS(dat.sub.impute.knn.lst, file = paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat.sub.impute.knn.lst.", Sys.Date(), ".rds"))
# jcheck <- dat.sub.impute.knn.lst$k4me1

dat.sub.impute.knn.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat.sub.impute.knn.lst.2022-05-06.rds")

cell2ptime.lsk <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/cell2ptime_lsk.2022-05-07.rds")
datlsklong.lsk <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat_lsk_long_ptime_lsk.2022-05-07.rds")

# Load PCA with arrows  ---------------------------------------------------


  # #CFB000
  # #00D1D0
  # #FF80F7

# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# cbPalette <- c("grey50", "#E69F00", "#0209b1", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette <- c("#CFB000", "#00D1D0", "#FC00F6")
jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew

eryth.traj <- c("HSCs", "ST", "LT", "MPPs", "MEP", "Eryths"); names(eryth.traj) <- eryth.traj
myeloid.traj <- c("HSCs", "ST", "LT", "MPPs", "CMP", "Basophils", "Granulocytes", "Monocytes", "DCs"); names(myeloid.traj) <- myeloid.traj
lymphoid.traj <- c("HSCs", "ST", "LT", "MPPs", "Bcells", "NKs"); names(lymphoid.traj) <- lymphoid.traj

traj.lst <- list("Eryths" = eryth.traj, 
                 "Granulocytes" = myeloid.traj, 
                 "Bcells" = lymphoid.traj)
jtypes <- names(traj.lst); names(jtypes) <- jtypes

jalphadots <- 1
jalphaarrows <- 0.5
jsizedots <- 1
jsizearrows <- 0.5
janglearrows <- 45
jlengtharrows <- 1.6

# get PC factors
pc1.factors.lst <- list("k27me3" = -1,
                        "k9me3" = -1, 
                        "k4me3" = 1,
                        "k4me1" = -1)

pc2.factors.lst <- list("k27me3" = 1,
                        "k9me3" = 1, 
                        "k4me3" = -1,
                        "k4me1" = 1)


# print(unique(dat.pca.merge.wide$ctype.from.LL))
ctypes.ordered <- c("Granulocytes", "Bcells", "Eryths", "NKs", "Monocytes", "DCs", "Basophils", "MEP", "CMP", "pDCs", "MPPs", "LT", "ST", "HSCs")
names(ctypes.ordered) <- ctypes.ordered

# Do K27me3 ---------------------------------------------------------------

# load PCA to get varexplained

jmark <- "k27me3"

for (jmark in jmarksnew){
  
  
  pc1.factor <- pc1.factors.lst[[jmark]]
  pc2.factor <- pc2.factors.lst[[jmark]]
  
  
  inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
  pca.out <- readRDS(inf.pca.out)
  inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
  load(inf.objs, v=T)
  
  dat.pca.merge.wide$ctype.from.LL <- factor(dat.pca.merge.wide$ctype.from.LL, levels = ctypes.ordered)
  
  varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
  varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)
  
  dat.arrows.pca <- GetArrows(dat.pca.merge.wide, grid.n = 30, jfactor = 1)
  m <- ggplot(mapping = aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, xend = pc1.factor * pc1.y, yend = pc2.factor * pc2.y)) +
    geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = jalphadots, size = jsizedots) +
    geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
    scale_color_identity() + 
    theme_bw() +
    ggtitle(jmark) + 
    xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
    ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  m <- ggplot(mapping = aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, xend = pc1.factor * pc1.y, yend = pc2.factor * pc2.y)) +
    geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcode), alpha = jalphadots, size = jsizedots) +
    # geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
    scale_color_identity() + 
    facet_wrap(~ctype.from.LL) + 
    theme_bw() +
    ggtitle(jmark) + 
    xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
    ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  # add staining 
  colors_and_legend <- tricolore::Tricolore(dat.sub.impute.knn.lst[[jmark]], 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
  dat.sub.impute.knn.lst[[jmark]]$rgb <- colors_and_legend$rgb
  
  
  dat.pca.merge.wide.annot <- left_join(dat.pca.merge.wide, subset(dat.sub.impute.knn.lst[[jmark]], select = c(cell, ckit_f, sca1_f, lin_f, rgb))) %>%
    rowwise() %>%
    mutate(rgb = ifelse(is.na(rgb), "grey85", rgb), 
           ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
    arrange(!is.na(sca1_f), sca1_f)
  
  # head(dat.pca.merge.wide.annot %>% arrange(desc(sca1_f)) %>% dplyr::select(cell, rgb, sca1_f, ckit_f, lin_f))
  # head(dat.pca.merge.wide.annot %>% arrange(desc(ckit_f)) %>% dplyr::select(cell, rgb, sca1_f, ckit_f, lin_f))
  # head(dat.pca.merge.wide.annot %>% arrange(desc(lin_f)) %>% dplyr::select(cell, rgb, sca1_f, ckit_f, lin_f))
  
  
  m <- ggplot(mapping = aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, xend = pc1.factor * pc1.y, yend = pc2.factor * pc2.y)) +
    geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
    geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
    scale_color_identity() + 
    theme_bw() +
    ggtitle(jmark) + 
    xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
    ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  m <- ggplot(mapping = aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, xend = pc1.factor * pc1.y, yend = pc2.factor * pc2.y)) +
    geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
    # geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
    scale_color_identity() + 
    theme_bw() +
    ggtitle(jmark) + 
    facet_wrap(~ctype.from.LL ) + 
    xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
    ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  
  dat.pca.merge.wide.annot.bytraj.lst <- lapply(jtypes, function(jtype){
    ctypes.keep <- traj.lst[[jtype]]
    subset(dat.pca.merge.wide.annot, ctype.from.LL %in% ctypes.keep) %>%
      mutate(traj = jtype)
  }) 
  
  for (jtype in jtypes){
    m <- ggplot(dat.pca.merge.wide.annot.bytraj.lst %>%
                  bind_rows() %>%
                  mutate(jalpha = ifelse(traj == jtype, 1, 0.25)), 
                aes(x = ptime.lsk, y = pc1.factor * pc1.x, color = colcode, alpha = jalpha, group = traj)) + 
      geom_point() + 
      ggtitle(paste(jmark, "Highlight", jtype)) + 
      scale_color_identity() + 
      scale_alpha_identity() + 
      geom_vline(xintercept = 0.5, linetype = "dotted") + 
      # facet_wrap(~ctype.from.LL) + 
      theme_bw() + 
      xlab("Sca1-cKit-Lin trajectory") + 
      theme(aspect.ratio=0.66, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(dat.pca.merge.wide.annot.bytraj.lst %>%
                  bind_rows() %>%
                  filter(!is.na(ptime.lsk)) %>%
                  mutate(jalpha = ifelse(traj == jtype, 1, 0.25)), 
                aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, color = colcode, alpha = jalpha, group = traj)) + 
      geom_point() + 
      ggtitle(paste(jmark, "Highlight", jtype)) + 
      scale_color_identity() + 
      scale_alpha_identity() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(dat.pca.merge.wide.annot.bytraj.lst %>%
                  bind_rows() %>%
                  filter(!is.na(ptime.lsk)) %>%
                  mutate(jalpha = ifelse(traj == jtype, 1, 0.25)), 
                aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, color = rgb, alpha = jalpha, group = traj)) + 
      geom_point() + 
      ggtitle(paste(jmark, "Highlight", jtype)) + 
      scale_color_identity() + 
      scale_alpha_identity() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    
    
    # # fit principal curve
    # library(princurve)
    # inmat <- subset(dat.pca.merge.wide.annot.bytraj.lst[["Eryths"]], 
    #                 select = c(cell, ptime.lsk, pc1.x, ctype.from.LL)) %>%
    #   dplyr::mutate(pc1.pretty = pc1.factor * pc1.x) %>%
    #   filter(!is.na(ptime.lsk)) %>%
    #   dplyr::select(-pc1.x) %>%
    #   filter(! (ctype.from.LL != jtype & ptime.lsk > 0.5))
    # 
    # ggplot(inmat, aes(x = ptime.lsk, y = pc1.pretty)) + 
    #   geom_point() + 
    #   geom_smooth(se = FALSE) + 
    #   theme_bw() + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    # inmat <- as.data.frame(inmat)
    # rownames(inmat) <- inmat$cell
    # inmat$cell <- NULL
    # inmat$ctype.from.LL <- NULL
    # inmat <- as.matrix(inmat)
    # 
    # 
    # 
    # plot(inmat)
    # pcurve.out <- principal_curve(x = inmat)
    # plot(pcurve.out$s)
    # points(inmat)
    # 
    }
  
  # show LSK
  jlong <- subset(datlsklong.lsk, mark == jmark, select = c(cell, sca1_f_impute, ckit_f_impute, lin_f_impute, ptime.lsk)) %>%
    data.table::melt(id.vars = c("cell", "ptime.lsk"), measure.vars = c("sca1_f_impute", "ckit_f_impute", "lin_f_impute"))
  
  m <- ggplot(jlong, aes(x = ptime.lsk, y = value, color = variable, group = variable)) + 
    geom_point(alpha = 0.05) + 
    geom_smooth(method = "gam", se = FALSE) + 
    theme_bw() + 
    xlab("Sca1-cKit-Lin trajectory") + 
    ylab("Relative marker levels") + 
    geom_vline(xintercept = 0.5, linetype = "dotted") + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  m <- ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = pc1.factor * pc1.x, color = colcode)) + 
    scale_color_identity() + 
    geom_point() + 
    geom_vline(xintercept = 0.5, linetype = "dotted") + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  for (jtype in jtypes){
    m <- ggplot(dat.pca.merge.wide.annot.bytraj.lst %>%
                  bind_rows() %>%
                  mutate(jalpha = ifelse(traj == jtype, 1, 0.2)) %>%
                  arrange(jalpha), 
                aes(x = ptime.lsk, y = pc2.factor * pc2.x, color = colcode, alpha = jalpha, group = traj)) + 
      geom_point() + 
      ggtitle(paste(jmark, "Highlight", jtype)) + 
      scale_color_identity() + 
      scale_alpha_identity() + 
      theme_bw() + 
      geom_vline(xintercept = 0.5, linetype = "dotted") + 
      xlab("Sca1-cKit-Lin trajectory") + 
      theme(aspect.ratio=0.66, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  
  
  m <- ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = pc2.factor * pc2.x, color = colcode)) + 
    scale_color_identity() + 
    geom_point() + 
    ggtitle(jmark) + 
    geom_vline(xintercept = 0.5, linetype = "dotted") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}
dev.off()


