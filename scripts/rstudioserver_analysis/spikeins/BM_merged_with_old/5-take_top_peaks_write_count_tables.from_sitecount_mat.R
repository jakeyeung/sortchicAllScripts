# Jake Yeung
# Date of Creation: 2020-11-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/5-take_top_peaks_write_count_tables.from_sitecount_mat.R
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
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDA from peaks  ----------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat"
assertthat::assert_that(dir.exists(indir))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMAllMerged2.from_peaks.sitecount_mat/filtNAcells_topbins"
dir.create(outdir)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

outlst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.K-30.Robj")
  inf <- file.path(indir, fname)
  print(inf)
  load(inf, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)
  return(list(tm.result = tm.result, dat.umap = dat.umap, count.mat = count.mat))
})


inf.rdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/LDA_downstream.BMround2all.merge_with_old.cleaner/LDA_downstream_objects.2020-11-05.again.RData"
load(inf.rdata, v=T)

tm.result.lst <- lapply(jmarks, function(jmark){
  tm.result <- outlst[[jmark]]$tm.result
})

count.mat.lst <- lapply(jmarks, function(jmark){
  count.mat <- outlst[[jmark]]$count.mat
})

dat.umap.lst <- lapply(jmarks, function(jmark){
  dat.umap <- outlst[[jmark]]$dat.umap %>%
    left_join(., subset(dat.merge, select = c(cell, cluster)))
  dat.umap$mark <- jmark
  return(dat.umap)
}) %>%
  bind_rows()


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.umap.lst %>% filter(mark == "H3K4me3"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")

ggplot(dat.umap.lst %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")



# Remove na cells ---------------------------------------------------------


cells.keep.lst <- lapply(jmarks, function(jmark){
  dat.tmp <- subset(dat.merge, mark == jmark)
  return(subset(dat.tmp, !is.na(cluster))$cell)
})
# cells.keep <- subset(dat.merge, !is.na(cluster))$cell


# Take topn ---------------------------------------------------------------

keeptop <- 250
bins.keep.lst <- lapply(jmarks, function(jmark){
  jtmp <- tm.result.lst[[jmark]]$terms
  jtopics <- rownames(jtmp)
  names(jtopics) <- jtopics
  topbins.lst <- lapply(jtopics, function(jtopic){
    xvec <- sort(tm.result.lst[[jmark]]$terms[jtopic, ], decreasing = TRUE)
    return(names(xvec)[1:keeptop])
  })
  topbins <- unique(unlist(topbins.lst))
  return(topbins)
})


# count tables ------------------------------------------------------------

count.mat.filt.lst <- lapply(jmarks, function(jmark){
  count.mat <- count.mat.lst[[jmark]]
  cells.keep.tmp <- colnames(count.mat) %in% cells.keep.lst[[jmark]]
  rows.keep.tmp <- rownames(count.mat) %in% bins.keep.lst[[jmark]]
  count.mat.filt <- count.mat[rows.keep.tmp, cells.keep.tmp]
  print(dim(count.mat))
  print(dim(count.mat.filt))
  return(count.mat.filt)
})


# Write -------------------------------------------------------------------


lapply(jmarks, function(jmark){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_topbins.rds"))
  saveRDS(count.mat.filt.lst[[jmark]], file = outf)
})

