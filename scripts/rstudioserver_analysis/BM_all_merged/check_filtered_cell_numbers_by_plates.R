# Jake Yeung
# Date of Creation: 2020-02-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/check_filtered_cell_numbers_by_plates.R
# Check 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)


jmark <- "H3K4me3"
inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-04.var_filt.UnenrichedAndAllMerged/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-04.Unenriched.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-04.Unenriched.K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)


# Show UMAP and count numbers ---------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


topics.mat <- posterior(out.lda)$topics

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

table(dat.umap.long$plate)

# Plot raw variance to check ----------------------------------------------

dat.var.raw <- CalculateVarRaw(as.matrix(count.mat), merge.size = 1000, calculate.ncuts = TRUE)

dat.merge <- left_join(dat.umap.long, dat.var.raw) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

ggplot(dat.merge, aes(x = ncuts, y = ncuts.var)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + facet_wrap(~plate)



