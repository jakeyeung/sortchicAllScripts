# Jake Yeung
# Date of Creation: 2020-03-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/Blood_with_intestines/explore_LDA_output_blood.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)


library(scchicFuncs)

library(hash)

library(hash)
library(igraph)
library(umap)

# Load quickly  -----------------------------------------------------------

# jmark <- "k27me3"
# jmark <- "k4me3"
jmark <- "k4me1"
# jmark <- "k9me3"m
inf1 <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestinesWithBlood.2020-02-29/lda_outputs.count_mat.blood_", jmark, ".countcutoff_1000-1000-1000-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.blood_", jmark, ".countcutoff_1000-1000-1000-1000-1000.TAcutoff_0.5.K-30.Robj")
# inf1 <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestinesWithBlood.2020-02-29/lda_outputs.count_mat.bloodIntest_k9me3.countcutoff_1000-1000-1000-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.bloodIntest_k9me3.countcutoff_1000-1000-1000-1000-1000.TAcutoff_0.5.K-30.Robj")
load(inf1, v=T)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
tm.result <- posterior(out.lda)
dat.umap.out <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))

ggplot(dat.umap.out, aes(x = umap1, y = umap2)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap.out, dat.var)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1)

# mark which is blood and intestine
dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(experi = ClipLast(as.character(cell), jsep = "_"))


ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)
