# Jake Yeung
# Date of Creation: 2022-05-17
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/34-check_lsk_on_PCA_fix_colors_ternary.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew

jalphadots <- 1
jalphaarrows <- 0.5
jsizedots <- 1
jsizearrows <- 0.5
janglearrows <- 45
jlengtharrows <- 1.6

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed"

outpdf <- file.path(outdir, paste0("LSK_metadata_colorcode_more_contrast_plots.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# Get new colors for celltypes ---------------------------------------------------------

dat.colors <- fread("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt")


# Load  -------------------------------------------------------------------



cell2ptime.lsk <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/cell2ptime_lsk.2022-05-07.rds")

dat.sub.impute.knn.lst <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/dat.sub.impute.knn.lst.2022-05-06.rds")

dat.sub.impute.knn.lst <- lapply(dat.sub.impute.knn.lst, function(jdat){
  jdat <- left_join(jdat, dat.colors)
})

ctypes.ordered <- c("Granulocytes", "Bcells", "Eryths", "NKs", "Monocytes", "DCs", "Basophils", "MEP", "CMP", "pDCs", "MPPs", "LT", "ST", "HSCs")
names(ctypes.ordered) <- ctypes.ordered

jmark <- "k27me3"

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

colors_and_legend <- tricolore::Tricolore(dat.sub.impute.knn.lst %>% bind_rows(), 
                                          p1 = 'sca1_f_impute', p2 = 'ckit_f_impute', p3 = 'lin_f_impute', 
                                          hue = 0.1, chroma = 1, lightness = 1, spread = 2)
print(colors_and_legend)

dat.sub.impute.knn.lst[[jmark]]$rgb <- colors_and_legend$rgb

dat.pca.merge.wide.annot <- left_join(dat.pca.merge.wide, subset(dat.sub.impute.knn.lst[[jmark]], select = c(cell, ckit_f, sca1_f, lin_f, rgb))) %>%
  rowwise() %>%
  mutate(rgb = ifelse(is.na(rgb), "grey85", rgb), 
         ptime.lsk = AssignHash(x = cell, jhash = cell2ptime.lsk, null.fill = NA)) %>%
  arrange(!is.na(sca1_f), sca1_f)




# Show UMAP  --------------------------------------------------------------

ggplot(dat.pca.merge.wide.annot, aes(x = umap1, y = umap2, color = rgb)) + 
  geom_point() + 
  scale_color_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(mapping = aes(x = 1 * pc1.x, y = -1 * pc2.x, xend = 1 * pc1.y, yend = -1 * pc2.y)) +
  geom_point(data = dat.pca.merge.wide.annot, mapping = aes(color = rgb), alpha = jalphadots, size = jsizedots) +
  geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows), data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
  scale_color_identity() + 
  theme_bw() +
  ggtitle(jmark) + 
  xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
  ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# Compare LSK ptime vs chrom ptime ----------------------------------------

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = -1 * pc1.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pca.merge.wide.annot, aes(x = ptime.lsk, y = pc2.x, color = colcode)) + 
  scale_color_identity() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Fix colors  -------------------------------------------------------------

dat.sub.impute.knn.rgb.lst <- lapply(jmarksnew, function(jmark){
  jdat <- dat.sub.impute.knn.lst[[jmark]] 
  colors_and_legend <- tricolore::Tricolore(jdat, 
                                            p1 = 'sca1_f_impute', p2 = 'ckit_f_impute', p3 = 'lin_f_impute', 
                                            hue = 0.1, chroma = 1, lightness = 1, spread = 2)
  jdat$rgb <- colors_and_legend$rgb
  jdat <- jdat %>%
    arrange(desc(-1 * sca1_f_impute + 1 * lin_f_impute))
  return(jdat)
})

m.lst <- lapply(jmarksnew, function(jmark){
  m <- ggplot(dat.sub.impute.knn.rgb.lst[[jmark]], aes(x = umap1, y = umap2, color = rgb)) + 
    geom_point() + 
    scale_color_identity() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst)

m.lst <- lapply(jmarksnew, function(jmark){
  m <- ggplot(dat.sub.impute.knn.rgb.lst[[jmark]], aes(x = umap1, y = umap2, color = colcodenew)) + 
    geom_point() + 
    scale_color_identity() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst)

dev.off()

for (jmark in jmarksnew){
  outftmp <- file.path(outdir, paste0("LSK_metadata_colorcode_more_contrast.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(x = dat.sub.impute.knn.rgb.lst[[jmark]], file = outftmp, sep = "\t")
}
saveRDS(colors_and_legend, file = file.path(outdir, paste0("LSK_color_and_legend_object.", Sys.Date(), ".rds")))

