# Jake Yeung
# Date of Creation: 2020-11-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_celltypes.R
# Get bone marrow celltypes

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDA output ---------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04"


outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/BM_UMAP_celltypes.", Sys.Date(), ".pdf")

pdf(file = outpdf, useDingbats = FALSE)

out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj")
  inf <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  tm.result <- posterior(out.lda)
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  return(list(dat.umap = dat.umap, tm.result = tm.result, count.mat = count.mat))
})


dat.umap.lda <- lapply(jmarks, function(jmark){
  jout <- out.lst[[jmark]]
  jdat <- jout$dat.umap
  jdat$mark <- jmark
  return(jdat)
}) %>%
  bind_rows()

ggplot(dat.umap.lda %>% filter(mark == "H3K4me1"), aes(x = umap1, y = umap2)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Load cell cluster -------------------------------------------------------

inf.rdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/LDA_downstream.BMround2all.merge_with_old.cleaner/LDA_downstream_objects.2020-11-05.again.RData"
load(inf.rdata, v=T)

dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "Eryth", "Eryths", cluster))
# merge
dat.umap.lda.annot <- left_join(dat.umap.lda, subset(dat.merge, select = c(cell, cluster), by = c("cell")))

ctypes <- sort(unique(subset(dat.umap.lda.annot, mark == "H3K4me1")$cluster))
ctypes2 <- gsub(pattern = "Bcells", replacement = "Lymphoid", ctypes)
# ctypes2 <- sort(unique(subset(dat.umap.lda.annot, mark == "H3K9me3")$cluster))

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jsub <- dat.umap.lda.annot %>% filter(mark == jmark) 
  if (jmark != "H3K9me3"){
    jsub$cluster <- factor(jsub$cluster, levels = ctypes)
    jctypes <- ctypes
  } else {
    jsub$cluster <- factor(jsub$cluster, levels = ctypes2)
    jctypes <- ctypes2
  }
  jpal <- cbPalette[sort(unique(as.numeric(jsub$cluster)))]
  jsub$jcol <- cbPalette[as.numeric(jsub$cluster)]
  m <- ggplot(jsub, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_minimal() + 
    scale_color_manual(values = jpal) + 
    # scale_color_identity(labels = unique(jsub$cluster), guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

print(m.lst)

dev.off()

