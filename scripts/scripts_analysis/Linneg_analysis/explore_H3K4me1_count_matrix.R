# Jake Yeung
# Date of Creation: 2019-08-15
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/explore_H3K4me1_count_matrix.R
# Explore count matrix of the new antibody 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(scchicFuncs)
library(irlba)
library(Matrix)
library(umap)
library(hash)
library(igraph)

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-Linneg/PZ-Bl6-BM-Linneg_H3K4me1_binfilt_cellfilt.2019-06-17.merge_with_B6.RData"
# inf2 <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_H3K4me1_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData"

load(inf, v=T)
count.h3k4me1.new <- count.dat$counts
# load(inf2, v=T)
# count.h3k4me1.old <- count.dat$counts

# add chr to old cnames
# rownames(count.h3k4me1.old) <- paste("chr", rownames(count.h3k4me1.old), sep = "")

# colnames(count.h3k4me1.new) <- paste("new_", colnames(count.h3k4me1.new), sep = "")
# colnames(count.h3k4me1.old) <- paste("old_", colnames(count.h3k4me1.old), sep = "")

# common.rows <- intersect(rownames(count.h3k4me1.new), rownames(count.h3k4me1.old))
# assertthat::assert_that(length(common.rows) > 0)
# count.h3k4me1.merged <- Matrix::cbind2(count.h3k4me1.new[common.rows, ], count.h3k4me1.old[common.rows, ])
count.h3k4me1.merged <- count.h3k4me1.new

linneg.cells <- grepl(pattern = "^PZ", colnames(count.h3k4me1.new))
count.h3k4me1.merged <- count.h3k4me1.new[, linneg.cells]

# Do SVD ------------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

lsi.out <- RunLSI(as.matrix(count.h3k4me1.merged), n.components = 50)

umap.out <- umap(lsi.out$u, config = jsettings)

dat.umap <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(batch = strsplit(cell, "-")[[1]][[1]])

ggplot(dat.umap, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# dat.louvain <- DoLouvain(lsi.out$u, custom.settings.louv = jsettings, dat.umap)

# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
# m <- PlotXYWithColor(dat.louvain, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette)
# print(m)



# label by something


# Find cell types  --------------------------------------------------------




