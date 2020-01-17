# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/fit_GLM_test.R
# Try fit GLM

rm(list=ls())

library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

# Load and fit ------------------------------------------------------------

inf.glm <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_test_regressvar_glmout.rds"
inf.glm2 <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_test_noregressvar_glmout.rds"
glm.out <- readRDS(inf.glm)
glm.out2 <- readRDS(inf.glm2)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

topics.mat <- glm.out$factors
topics.mat.noregress <- glm.out2$factors

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out <- umap(topics.mat.noregress, config = jsettings)
dat.umap.long.noregress <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long.noregress <- DoLouvain(topics.mat.noregress, jsettings, dat.umap.long.noregress)
ggplot(dat.umap.long.noregress, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate()

# load old louvain

inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
load(inf, v=T)

tm.result <- posterior(out.lda)
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

topics.mat.orig <- tm.result$topics
umap.out <- umap(topics.mat.orig, config = jsettings)
dat.umap.long.orig <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)

dat.umap.long.orig <- DoLouvain(topics.mat.orig, jsettings, dat.umap.long.orig, clstr.cname = "old_louvain")

dat.umap.long.orig <- left_join(dat.umap.long.orig, dat.var)


dat.merge <- left_join(dat.umap.long, dat.var)
dat.merge.noregress <- left_join(dat.umap.long.noregress, dat.var)
dat.merge <- left_join(dat.merge, subset(dat.umap.long.orig, select = c(cell, old_louvain)))
dat.merge.noregress <- left_join(dat.merge.noregress, subset(dat.umap.long.orig, select = c(cell, old_louvain)))

m.old <- ggplot(dat.umap.long.orig, aes(x = umap1, y = umap2, color = old_louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + ggtitle("k4me1 orig")
m.new <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = old_louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + ggtitle("k4me1 after")

JFuncs::multiplot(m.old, m.new)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle("k4me1 after var correction")

ggplot(dat.merge.noregress, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.umap.long.orig, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# check plates
dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(cell, jsep = "-")) 

dat.merge.noregress <- dat.merge.noregress %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(cell, jsep = "-")) 

dat.umap.long.orig <- dat.umap.long.orig %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(cell, jsep = "-"))

ggplot(dat.merge, aes(x = umap1, y = umap2, color = old_louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)  + 
  facet_wrap(~plate) + ggtitle("Corrected for variance k4me1")

ggplot(dat.merge.noregress, aes(x = umap1, y = umap2, color = old_louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)  + 
  facet_wrap(~plate)

ggplot(dat.umap.long.orig, aes(x = umap1, y = umap2, color = old_louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~plate) + ggtitle("Orig k4me1 includes different antibodies, bad linneg, everything")



# check X fit with variance? 
# X <- matrix(dat.var$cell.var.within.sum.norm, ncol = 1, byrow = TRUE)
# dat.X <- as.matrix(glm.out$coefX)[, 2] %*% t(X)

