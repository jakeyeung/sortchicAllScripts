# Jake Yeung
# Date of Creation: 2020-08-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/celltyping_with_chix/2-LDA_projections_downstream.R
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
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "H3K27me3"


# Load annots for double chic ---------------------------------------------



# Load annots -------------------------------------------------------------

inf.annot <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables/BM_AllMerged.", jmark, ".cell_cluster_table.txt"))
dat.annot <- fread(inf.annot)

# Load LDA projections ----------------------------------------------------

inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_projections/project_unmixed_K27m3.project_on_AllMerged_", jmark, ".RData"))
load(inf.lda, v=T)

topics.mat <- posterior(out.objs$out.lda)$topics
topics.mat.proj <- out.lda.predict$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.orig <- umap(topics.mat, jsettings)
umap.proj <- predict(umap.orig, data = topics.mat.proj)

dat.umap.proj <- data.frame(cell = rownames(umap.proj), umap1 = umap.proj[, 1], umap2 = umap.proj[, 2], type = "proj", stringsAsFactors = FALSE)

# umap.proj.alone <- predict(umap.orig, data = topics.mat.proj)

dat.umap.orig <- DoUmapAndLouvain(topics.mat, jsettings = jsettings)

dat.umap.orig$type <- "orig"


ggplot(dat.umap.orig, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.orig.annot <- dat.umap.orig %>%
  rowwise() %>%
  mutate(mark = GetMarkFromStr(cell),
         cond = GetCondFromSamp(samp = cell, mark = mark))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.orig.annot, aes(x = umap1, y = umap2, color = cond)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.merged <- rbind(subset(dat.umap.orig, select = c(cell, umap1, umap2, type)), subset(dat.umap.proj, select = c(cell, umap1, umap2, type)))

ggplot(dat.umap.merged, aes(x = umap1, y = umap2, color = type)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Mix orig and proj count mats and rerun LDA  -----------------------------









