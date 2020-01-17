# Jake Yeung
# Date of Creation: 2019-12-20
# File: ~/projects/scchic/scripts/scripts_analysis/intestinal_analysis/explore_old_LDA_analysis.R
# Check old LDA

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

library(topicmodels)
library(irlba)

# Load data ---------------------------------------------------------------

jmark <- "k9me3"

jmark <- "k4me3"
indir <- "/home/jyeung/hpc/intestinal_scchic/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_intestinal_OUD3907/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.mapq_60"
fname <- paste0("lda_out_meanfilt.intestinal_OUD3907_", jmark, "_binfilt_cellfilt.CountThres0.K-20.Robj")
inf <- file.path(indir, fname)

assertthat::assert_that(file.exists(inf))

load(inf, v=T)



# Get umap  ---------------------------------------------------------------

tm.result <- posterior(out.lda[[1]])
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(plate = strsplit(cell, "-")[[1]][[6]])
print(head(dat.umap.long))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = plate)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle(jmark)

# check raw mat
lsi.out <- RunLSI(as.matrix(count.mat))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out.lsi <- umap(lsi.out$u, config = jsettings)
dat.umap.long.lsi <- data.frame(cell = rownames(umap.out.lsi$layout), umap1 = umap.out.lsi$layout[, 1], umap2 = umap.out.lsi$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(plate = strsplit(cell, "-")[[1]][[6]])

ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2, color = plate)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle(jmark)


nnzero(count.mat) / length(count.mat)

