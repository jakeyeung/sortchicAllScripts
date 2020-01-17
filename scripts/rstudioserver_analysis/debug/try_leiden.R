# Jake Yeung
# Date of Creation: 2020-01-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/try_leiden.R
# Try leiden

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


library(leidenbase)
library(jclustering)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)
library(JFuncs)

DoLouvain()

DoLeiden <- function(topics.mat, K=30, res_param=10^seq(-5, 0), dat.umap.long = NULL, random_seed = 123, weight=FALSE, verbose = TRUE){
  pdat <- data.frame(cell = rownames(topics.mat))
  rownames(pdat) <- pdat$cell
  rownames(pdat) <- pdat$cell
  l.out <- leiden_clustering(topics.mat, pdat, random_seed = random_seed, weight = weight, verbose=verbose, resolution_parameter = res_param, k = K)
  dat.leiden <- data.frame(cell = names(l.out$optim_res$membership), cluster = paste("leiden", as.character(l.out$optim_res$membership), sep = "_"))
  if (!is.null(dat.umap.long)){
    dat.leiden <- left_join(dat.umap.long, dat.leiden)
  }
  return(dat.leiden)
}


topics.mat <- readRDS("~/data/topics_mat.rds")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


data <- as.data.frame(topics.mat)
# colnames(data) <- paste("gene", seq(ncol(data)), sep = "_")
pd <- data.frame(cell = rownames(data))
rownames(pd) <- pd$cell

l.out <- leiden_clustering(data, pd, random_seed = 123, weight = FALSE, verbose=TRUE, resolution_parameter = c(0.0001, 0.001, 0.01, 0.1, 1), k = jsettings$n_neighbors)

dat.leiden <- data.frame(cell = names(l.out$optim_res$membership), leiden = as.character(l.out$optim_res$membership))

dat.merge <- left_join(dat.umap.long, dat.leiden)

m.louv <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
m.leid <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = leiden)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

multiplot(m.louv, m.leid, cols = 1)

