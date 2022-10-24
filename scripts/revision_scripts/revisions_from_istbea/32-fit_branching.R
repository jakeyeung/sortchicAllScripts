# Jake Yeung
# Date of Creation: 2022-05-06
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/32-fit_branching.R
# 

rm(list=ls())

library(igraph)
library(mgcv)
library(quadprog)
library(pcaMethods)
library(Rcpp)
library(inline)
library(RcppArmadillo)
library(pbapply)
library(crestree)
library(ggplot2); library(gridExtra); library(grid)

jmark <- "k4me3"
inf.rds <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits/signal_extrap_vectors_regions_ctypes_from_args.", jmark, ".rds")
dat.rds <- readRDS(inf.rds)

mat.signal <- lapply(dat.rds, function(jdat){
  jdat$signal.vec
}) 
mat.signal <- do.call(rbind, mat.signal)

mat.input <- mat.signal

mat.weights <- matrix(1, nrow = nrow(mat.input), ncol = ncol(mat.input))
rownames(mat.weights) <- rownames(mat.input)
colnames(mat.weights) <- colnames(mat.input)
mat.weights <- sweep(mat.weights, MARGIN = 2, STATS = colSums(mat.weights), FUN = "/")

M <- nrow(mat.input)
lambda <- 250;
sigma <- 0.04


inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
pca.out <- readRDS(inf.pca.out)
inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
load(inf.objs, v=T)

emb <- as.data.frame(subset(dat.pca.merge.wide, select = c(cell, pc1.x, pc2.x)))
rownames(emb) <- dat.pca.merge.wide$cell
emb$cell <- NULL
emb <- as.matrix(emb)
  
system.time(
  ppt <- ppt.tree(X=mat.input, W=mat.weights, emb=emb, lambda=250, sigma=0.04, metrics="cosine", M=M, err.cut = 5e-3, n.steps=30, seed=1, plot=FALSE)
)

plotppt(ppt,emb,tips=FALSE,cex.tree = 0.1,cex.main=0.2,lwd.tree = 1)
plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)







# Check -------------------------------------------------------------------



data(crest)
emb <- crest$emb
clcol <- crest$clcol
nc.cells <- crest$nc.cells
wgm <- crest$wgm
wgwm <- crest$wgwm # matrix of expression weights
# fpm <- read.table("http://pklab.med.harvard.edu/ruslan/neural_crest/fpm.txt",header=TRUE)
# fpm <- as.matrix(fpm)
# write.table(fpm, file = "/nfs/scistore12/hpcgrp/jyeung/IST_data/public_data/neural_crest/neural_crest_Soldatov_2019_fpm.txt",
#             quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
fpm <- read.table("/nfs/scistore12/hpcgrp/jyeung/IST_data/public_data/neural_crest/neural_crest_Soldatov_2019_fpm.txt.gz")
# save table locally
fpm <- as.matrix(fpm)
genes.tree <- crest$genes.tree







# Fit tree ----------------------------------------------------------------

M <- length(nc.cells);
lambda <- 250;
sigma <- 0.04

system.time(
  ppt <- ppt.tree(X=wgm[,nc.cells], W=wgwm[,nc.cells], emb=emb, lambda=250, sigma=0.04, metrics="cosine", M=M, err.cut = 5e-3, n.steps=30, seed=1, plot=FALSE)
)

plotppt(ppt,emb,tips=FALSE,cex.tree = 0.1,cex.main=0.2,lwd.tree = 1)
plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)

# have to do some manual cleanups iteratively and then assign root
# ppt <- cleanup.branches(ppt,tips.remove = c(222))
ppt <- setroot(ppt, 188)
ppt <- project.cells.onto.ppt(ppt,n.mapping = 100)


