# Jake Yeung
# Date of Creation: 2022-01-18
# File: ~/projects/scchic/scripts/macbook_analysis_2021/new_experiments/1-check_LDA_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)


library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

# Load --------------------------------------------------------------------


jmarks <- c("k27me3"); names(jmarks) <- jmarks

# jmark <- "k4me3"
# jmark <- "k9me3"
# inf <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/new_experiments/LDA_outputs/K562/ldaAnalysis_fripfilt_K562_", jmark, "/lda_outputs.count_mat_merged_with_old.K562_", jmark, ".0.8_0.5_3000.2022-01-10/ldaOut.count_mat_merged_with_old.K562_", jmark, ".0.8_0.5_3000.2022-01-10.Robj")

# jmark <- "k27me3"

outpdf <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_macbook_2021/BM_new_experiment/BM_new_experiments.merge_with_old.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)
for (jmark in jmarks){
  
  # jmark <- "k4me1-k9me3"
  inf <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/new_experiments/LDA_outputs/ldaAnalysis_fripfilt_BM_", jmark, "/lda_outputs.count_mat_merged_with_old.", jmark, ".2022-01-09/ldaOut.count_mat_merged_with_old.", jmark, ".2022-01-09.Robj")
  # load(inf, v=T)
  load(inf, v=T)
  
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  
  dat.umap <- DoUmapAndLouvain(topics.mat = tm.result$topics, jsettings = jsettings)
  
  dat.umap.annot <- dat.umap %>%
    rowwise() %>%
    mutate(plate = ClipLast(x = cell, jsep = "_"),
           experi = strsplit(cell, split = "-")[[1]][[2]])
  
  # Annotate old and new  ---------------------------------------------------
  
  m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point() + 
    ggtitle(jmark) + 
    facet_wrap(~experi) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}
dev.off()



