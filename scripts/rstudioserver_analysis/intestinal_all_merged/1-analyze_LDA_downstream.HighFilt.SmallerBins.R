# Jake Yeung
# Date of Creation: 2020-01-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-analyze_LDA_downstream.HighFilt.SmallerBins.R
# High Filt and Smaller Bins


library(scchicFuncs)
library(JFuncs)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(igraph)
library(umap)
library(topicmodels)



# Load LDA  ---------------------------------------------------------------

jwinbase <- "10000_5000"
# jwin <- "20000_10000"
jwin <- "50000_25000"
jmark <- "k36me3"
jprefix <- "AllMerged"
inf <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2020-01-11.", 
              jwinbase, "/lda_outputs.mat.Scraped.", jprefix, ".", jmark, 
              ".TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_", jwin, 
              ".2020-01-11.K-30.binarize.FALSE/ldaOut.mat.Scraped.", jprefix, ".", jmark, 
              ".TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_", jwin, ".2020-01-11.K-30.Robj")

inf.check <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2020-01-11.10000_5000/lda_outputs.mat.Scraped.AllMerged.k4me3.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_20000_10000.2020-01-11.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.k4me3.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_20000_10000.2020-01-11.K-30.Robj"
assertthat::assert_that(file.exists(inf))


load(inf, v=T)

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


