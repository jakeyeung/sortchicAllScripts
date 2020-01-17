# Jake Yeung
# Date of Creation: 2020-01-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/show_var_across_systems_use_objs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(parallel)
library(hash)
library(igraph)
library(umap)


# Load data ---------------------------------------------------------------

inf <- "/home/jyeung/data/from_rstudioserver/intestinalchic/var_across_intestinal_marks.rds"

out.lst <- readRDS(inf)

# out.lst$k9me3 <- NULL
# get fraction of reads in chr19 and chr11?

chrs.high <- c("chr11", "chr19")
frac.reads <- lapply(out.lst, function(x){
  jdat <- x$reads.by.chromo.sum %>%
    rowwise() %>%
    mutate(chromo.tmp = chromo %in% chrs.high) %>%
    group_by(cell, chromo.tmp) %>%
    summarise(ncuts = sum(ncuts)) %>%
    group_by(cell) %>%
    mutate(nfrac = ncuts / sum(ncuts))
  return(jdat)
}) %>% 
  bind_rows()

frac.reads.filt <- subset(frac.reads, chromo.tmp)

ggplot(out.lst$k27me3$reads.by.chromo.sum, aes(x = ncuts)) + geom_density() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_log10() + facet_wrap(~chromo)

# remove heterochromatin

# Combine counts across chromos and then do PCA  --------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
chromos.keep <- jchromos[!jchromos %in% c("chrX", "chrY")]

mats.combined <- lapply(out.lst, function(x){
  mat <- tidyr::spread(data = x$reads.by.chromo %>% select(-label), cell, ncuts)
  rownames(mat) <- mat$chromo
  mat <- subset(mat, chromo %in% chromos.keep)
  mat$chromo <- NULL
  # colnames(mat) <- paste(jlab, colnames(mat), sep = "_")
  return(as.matrix(mat))
}) 
mats.combined <- do.call(cbind, mats.combined)

# remove sex chromosomes

dat.var.combined <- lapply(out.lst, function(x){
  jlab <- unique(x$reads.by.chromo$label)
  assertthat::assert_that(length(jlab) == 1)
  x$dat.meta$mark <- jlab
  return(x$dat.meta)
}) %>%
  bind_rows()

mat.input <- scale(t(mats.combined), center = TRUE, scale = TRUE)

pca.out <- prcomp(mat.input, center = FALSE, scale. = FALSE)

dat.pca <- data.frame(cell = rownames(pca.out$x), pca.out$x, stringsAsFactors = FALSE) %>%
  left_join(., dat.var.combined) %>%
  left_join(., subset(frac.reads.filt, select = c(-chromo.tmp, -ncuts)))

ggplot(dat.pca, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~mark)

ggplot(dat.pca, aes(x = PC1, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~mark)

ggplot(dat.pca, aes(x = PC2, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~mark)

ggplot(dat.pca, aes(x = nfrac, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~mark)

# what is PC2?
dat.loadings <- data.frame(chromo = rownames(pca.out$rotation), pca.out$rotation, stringsAsFactors = FALSE)

# library(ggrepel)
ggplot(dat.loadings, aes(x = PC1, y = PC2, label = chromo)) + geom_point() + geom_text_repel() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())






