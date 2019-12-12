# Jake Yeung
# Date of Creation: 2019-12-08
# File: ~/projects/dblchic/scripts/macbook_analysiis/create_louvain_objects_from_LDA.R
# Create louvain objects from LDA, then fit multinomials 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

# Load LDA outputs --------------------------------------------------------


jwin <- "20000"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo.other_winsizes/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.winsize_bsize_", jwin, ".K-30.binarize.TRUE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.winsize_bsize_", jwin, ".K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell),
         experi = gsub("-G2$", "", experi))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~experi)



# Load raw  ---------------------------------------------------------------


