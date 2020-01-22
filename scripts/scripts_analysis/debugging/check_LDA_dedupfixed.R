# Jake Yeung
# Date of Creation: 2019-12-02
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_dedupfixed.R
# LDA dedup fixed

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)


# Load data  --------------------------------------------------------------

jmark <- "H3K4me3"
jsuff <- "UnenrichedXLinneg"
jbin <- "TRUE"

# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_", 
#               jsuff, "_", jmark, "_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.binarize.", jbin, 
#               "/ldaOut.B6BM_", jsuff, "_", jmark, "_pcutoff_0.5_binfilt_cellfilt.2019-11-30.K-30.Robj")

inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed/lda_outputs.B6BM_",  jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.binarize.", jbin, "/ldaOut.B6BM_", jsuff, "_", jmark, "_pcutoff_TAcutoff_0.5.cellsize_1000.oldbins_binfilt_cellfilt.2019-12-01.K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)  %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell))
  
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# do variance

dat.impute.log <- log2(t(topics.mat %*% terms.mat))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.long.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)

ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check whether we see anything in the raw data? 
cells.linneg <- grepl(pattern = "Linneg", colnames(count.mat.orig))
cells.wt <- grepl(pattern = "^B6-13W1", colnames(count.mat.orig))

count.mat.linneg <- count.mat.orig[, cells.linneg]
count.mat.wt <- count.mat.orig[, cells.wt]




# # use entropy as a measurement??
# 
# 

# # Merge by louvains and plot  ---------------------------------------------
# 
# dat.umap.long.merge.louv <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.merge)
# 
# # visualize the variance??
# nnzero(count.mat.linneg) / length(count.mat.linneg)
# nnzero(count.mat.sc) / length(count.mat.sc)



