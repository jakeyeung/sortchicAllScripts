# Jake Yeung
# Date of Creation: 2019-11-04
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA.R
# Zebrafish

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

# Load dat ----------------------------------------------------------------

# jmark <- "H3K4me1"
jmark <- "H3K4me3"
# jmark <- "H3K9me3"
jbin <- "FALSE"

inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# Plot data ---------------------------------------------------------------

topics.mat <- posterior(out.lda[[1]])$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("CD41plus", cell))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)


# Do variance calculation -------------------------------------------------

jfac <- 10^6
jpseudo <- 0
dat.impute.log <- log2(t(posterior(out.lda[[1]])$topics %*% posterior(out.lda[[1]])$terms) * jfac + jpseudo)
jchromos.num <- seq(25)
jchromos <- paste("chr", jchromos.num, sep = "")

cells.var.chromo.merged <- CalculateVarAll(dat.impute.log, jchromos)


# Plot with variance  -----------------------------------------------------

dat.umap.long.merge <- left_join(dat.umap.long, cells.var.chromo.merged)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  scale_colour_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark) + 
  facet_wrap(~is.stem)

# do histogram?
ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm)) + geom_histogram() + facet_wrap(~is.stem, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark) + xlab("Intrachromosomal Variance")





