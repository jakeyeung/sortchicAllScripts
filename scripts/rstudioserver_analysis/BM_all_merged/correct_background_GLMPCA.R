# Jake Yeung
# Date of Creation: 2020-01-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/correct_background_GLMPCA.R
# Correct background GLMPCA 

rm(list=ls())


jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

# Load K4me3 --------------------------------------------------------------

jmark <- "H3K4me3"
inf <-  paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result <- posterior(out.lda)

outdir <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs"
assertthat::assert_that(dir.exists(outdir))
outf <- file.path(outdir, "fit_GLM_regress_intrachromvar_init_with_LDA.RData")
assertthat::assert_that(!file.exists(outf))

# Calculate intrachrom var  -----------------------------------------------

dat.var.raw <- CalculateVarRaw(as.matrix(count.mat), merge.size = 1000, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)

# calculate var from LDA

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

# Check dat.var.raw recapitulates dat.var ---------------------------------

dat.var.merge <- left_join(dat.var.raw, dat.var)

jcutoff <- 0.3  # ncuts.var 

ggplot(dat.var.merge, aes(y = cell.var.within.sum.norm, x = ncuts.var, color = cell.var.within.sum.norm)) + geom_point() +
  scale_x_log10() +  scale_y_log10() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + geom_vline(xintercept = jcutoff)


ggplot(dat.var.merge, aes(y = ncuts, x = ncuts.var, color = cell.var.within.sum.norm)) + geom_point() +
  scale_x_log10() +  scale_y_log10() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + geom_vline(xintercept = jcutoff)

ggplot(dat.var.merge, aes(y = ncuts, x = ncuts.var, color = cell.var.within.sum.norm)) + geom_point() +
  theme_bw() + scale_x_log10() + scale_y_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + geom_vline(xintercept = jcutoff)

# Refit GLM PCA  ----------------------------------------------------------

# filter cells below linear region? 
# cells.keep <- subset(dat.var.merge, ncuts.var > jcutoff)$cell

bins.keep <- 100

Y <- count.mat

assertthat::assert_that(all(rownames(Y) == colnames(tm.result$terms)))

# keep variable bins
bins.high.i <- as.data.frame(apply(tm.result$terms, MARGIN = 1, function(jcol) order(jcol, decreasing = TRUE)[1:bins.keep])) %>%
  unlist()

bins.high <- unique(rownames(Y)[bins.high.i])

Y.filt <- Y[bins.high, ]

var.vec <- dat.var.merge$ncuts.var
names(var.vec) <- dat.var.merge$cell

X <- data.frame(ncuts.var = var.vec, cell = names(var.vec))  # columns of 1s are implicit
X.reorder <- X[match(colnames(Y.filt), X$cell), ]
X.mat <- matrix(data = X.reorder$ncuts.var, ncol = 1, byrow = TRUE, dimnames = list(X.reorder$cell, "ncuts.var"))

# factors (U) is c by k
# loadings (V) is g by k
U.init <- scale(tm.result$topics, center = TRUE, scale = FALSE)  # c by k
V.init <- scale(tm.result$terms[, bins.high], center = TRUE, scale = FALSE)  # k by g
V.init <- t(V.init)  # now its g by k


assertthat::assert_that(all(rownames(X.mat) == colnames(Y.filt)))
assertthat::assert_that(all(rownames(U.init) == colnames(Y.filt)))
assertthat::assert_that(all(rownames(V.init) == rownames(Y.filt)))

# # check inits are sane
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# imput.mat <- V.init %*% t(U.init)
# pca.out <- prcomp(t(imput.mat), center = FALSE, scale. = FALSE)
# dat.pca <- data.frame(cell = rownames(pca.out$x), PC1 = pca.out$x[, 1], PC2 = pca.out$x[, 2], stringsAsFactors = FALSE)
# ggplot(dat.pca, mapping = aes(x = PC1, y = PC2)) + geom_point() 
# umap.out <- umap(pca.out$x, config = jsettings)
# dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long <- DoLouvain(pca.out$x, jsettings, dat.umap.long)
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
# dat.umap.long <- dat.umap.long %>%
#   rowwise() %>%
#   mutate(experi = ClipLast(cell, jsep = "_"))
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)



# # # test that multiplying the two vectors gives centered output?
# impute.gene.cell <- U.init[1, ] %*% V.init[, 1]
# impute.gene <- U.init %*% V.init[, 1]


system.time(
  glm.out <- glmpca(Y = Y.filt, L = 30, fam = "mult", ctl = list(maxIters = 250, eps = 1e-4), penalty = 1, verbose = FALSE, 
                    init = list(factors = U.init, loadings = V.init), X = X.mat, Z = NULL)
)
save(glm.out, Y, U.init, V.init, X.mat, out.lda, file = outf)
print(Sys.time() - jstart)




