# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/2-prepare_exprs_mat_from_GLMPCA.R
# Explore GLMPCA output and make exprs mat


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"



jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load annots -------------------------------------------------------------

dat.annot.lst <- lapply(jmarks, function(jmark){
  inf.annot <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt"))
  assertthat::assert_that(file.exists(inf.annot))
  dat.annot <- fread(inf.annot)
  dat.annot$mark <- jmark
  return(dat.annot)
})


# Load GLMPCA output ------------------------------------------------------


jmark.ref <- "H3K27me3"
# jmark.ref <- "H3K4me1"
# jmark.ref <- "H3K4me3"

# inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.remove_eryth/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.GLMPCA_var_correction.mergebinsize_1000.binskeep_250.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.5.winsorize_TRUE.2020-04-29.RData")
# 
# inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.remove_eryth/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.remove_eryth.GLMPCA_var_correction.mergebinsize_1000.binskeep_250.covar_cell.var.within.sum.norm.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-04-29.RData")
# 
# inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.GLMPCA_var_correction.mergebinsize_1000.binskeep_250.covar_cell.var.within.sum.norm.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-04-29.RData")
# 
# inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.poi.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.GLMPCA_var_correction.poi.mergebinsize_1000.binskeep_250.covar_cell.var.within.sum.norm.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-04-29.RData")

# inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.NoCorrection.poi.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.GLMPCA_var_correction.poi.mergebinsize_1000.binskeep_250.covar_cell.var.within.sum.norm.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-04-29.RData")
inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.Faster.poi.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.GLMPCA_var_correction.poi.mergebinsize_1000.binskeep_250.covar_cell.var.within.sum.norm.log2.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-04-29.RData")
# outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.Faster.poi.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks/ZF_", jmark.ref, ".AllMerged.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks.GLMPCA_var_correction.poi.mergebinsize_1000.binskeep_250.covar_cell.var.within.sum.norm.log2.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-04-29.mat")
# assertthat::assert_that(!file.exists(outf))


print(inf)

assertthat::assert_that(file.exists(inf))
load(inf, v=T)

# dat.var <- subset(dat.merge2, select = c(cell, cell.var.within.sum.norm.CenteredAndScaled, cell.var.within.sum.norm, cell.var.within.sum.norm.log2.CenteredAndScaled))
dat.var <- subset(dat.merge2, select = c(cell, cell.var.within.sum.norm.CenteredAndScaled, cell.var.within.sum.norm))

# Do uMAP from GLMPCA -----------------------------------------------------


dat.impute.glm <- as.matrix(glm.out$loadings) %*% as.matrix(t(glm.out$factors))

dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)

dat.merge <- left_join(dat.umap, dat.annot.lst[[jmark.ref]], by = "cell") %>%
  left_join(., dat.var, by = "cell")

ggplot(dat.merge, aes(x = umap1.x, y = umap2.x, color = louvain.x)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = umap1.x, y = umap2.x, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = umap1.x, y = umap2.x, color = log(cell.var.within.sum.norm))) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = umap1.x, y = umap2.x, color = log(cell.var.within.sum.norm))) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) +  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot an eryth gene


dat.glm <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
  left_join(., dat.annot.lst[[jmark.ref]], by = "cell") %>%
  left_join(., dat.var, by = "cell")

ggplot(dat.glm, aes(x = dim1, y = dim2, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.glm, aes(x = dim1, y = dim2, color = log(cell.var.within.sum.norm))) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.glm, aes(x = dim2, y = dim3, color = cluster)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

# top loadings in dim1
ldat <- data.frame(peaks = rownames(glm.out$loadings), loading = glm.out$loadings$dim1, stringsAsFactors = FALSE) %>%
  # arrange(desc(loading))
  arrange(loading)
(jpeak <- ldat$peaks[[1]])
impute.vec <- data.frame(cell = colnames(dat.impute.glm), imputed = dat.impute.glm[jpeak, ], stringsAsFactors = FALSE)


dat.glm2 <- left_join(dat.glm, impute.vec)

ggplot(dat.glm2, aes(x = dim1, y = dim2, color = imputed))  + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.merge, aes(x = umap1.y, y = umap2.y, color = cluster)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.glm2, aes(x = dim1, y = dim2, color = cluster))  + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark.ref)

ggplot(dat.glm2, aes(x = dim1, y = dim2, color = plate))  + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark.ref)

# ggplot(dat.glm2, aes(x = dim2, y = dim3, color = cluster))  + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_manual(values = cbPalette) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark.ref)


# Write inferred model to output ------------------------------------------

# exprs.mat <- data.frame(Gene.ID = rownames(dat.impute.glm), dat.impute.glm, stringsAsFactors = FALSE)
# 
# fwrite(x = exprs.mat, file = outf, sep = "\t")
# 
# plot(density(colSums(dat.impute.glm)))
# plot(density(rowSums(dat.impute.glm)))

# dat.impute.glm2 <- predict(object = glm.out)
# 
# plot(dat.impute.glm[, 1], dat.impute.glm2[, 1], pch = ".")
# 
# plot(dat.impute.glm[, 1], log(dat.impute.glm2[, 1]), pch = ".")
# 
# plot(exp(dat.impute.glm[, 1]), dat.impute.glm2[, 1], pch = ".")

# cmd <- paste0("gzip ")
# system()

