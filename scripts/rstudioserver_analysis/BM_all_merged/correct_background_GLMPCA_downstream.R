# Jake Yeung
# Date of Creation: 2020-01-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/correct_background_GLMPCA_downstream.R
# Check downstream

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.RData"
# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.RData"
# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.sizefactor.again.RData"
# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.UseOrigSZ_TRUE.niter_500.log_FALSE.svdinit_TRUE.RData"
# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.UseOrigSZ_TRUE.niter_500.log_FALSE.svdinit_TRUE.RData"
# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.UseOrigSZ_TRUE.niter_1000.log_FALSE.svdinit_TRUE.RData"

niter <- "250"
inf <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/fit_GLM_regress_intrachromvar_init_with_LDA.UseOrigSZ_TRUE.niter_", niter, ".log_FALSE.svdinit_TRUE.RData")
load(inf, v=T)


topics.mat <- glm.out$factors

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette <- rep(cbPalette, 2)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  ggtitle(paste("N iterations:", niter))

# recalculate intrachrom var???

dat.umap.long.merge <- left_join(dat.umap.long, dat.var.merge) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"))

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) +
# ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = ncuts.var)) +  
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  facet_wrap(~experi)
  
ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  ggtitle(paste("N iterations:", niter)) + facet_wrap(~experi)


ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm, y = ncuts.var)) +  
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


