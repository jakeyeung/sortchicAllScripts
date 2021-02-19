# Jake Yeung
# Date of Creation: 2021-01-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/17-LDA_checks_H3K4me1_H3K9me3_dbl.from_bins.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

jdate <- "2021-01-31"
# jmark <- "H3K4me1"
jsuffix <- "from_genomewide.RemoveHSPCs"
jsuffix2 <- "50kb_genomewide.RemoveHSPCs"

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks.withdbl <- c("H3K4me1", "H3K9me3", "H3K4me1-H3K9me3"); names(jmarks.withdbl) <- jmarks.withdbl


# outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/lda_and_clusters.", jsuffix2)
# dir.create(outdir)

# Metas  ------------------------------------------------------------------

indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})

# LDA  --------------------------------------------------------------------

indir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt.k4_k9_dynamic_bins_", jsuffix)
assertthat::assert_that(dir.exists(indir))


outs.lst <- lapply(jmarks.withdbl, function(jmark){
  print(jmark)
  dname <- paste0("lda_outputs.count_name.", jmark, ".k4_k9_", jsuffix2, ".", jdate, ".K-30.binarize.FALSE/ldaOut.count_name.", jmark, ".k4_k9_", jsuffix2, ".", jdate, ".K-30.Robj")
  inf <- file.path(indir, dname)
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

out.lda.lst <- lapply(jmarks.withdbl, function(jmark){
  jout.lda <- outs.lst[[jmark]]$out.lda
  return(jout.lda)
})

count.mat.lst <- lapply(jmarks.withdbl, function(jmark){
  jout.lda <- outs.lst[[jmark]]$count.mat
  return(jout.lda)
})

dat.umap.lst <- lapply(jmarks.withdbl, function(jmark){
  out.lda <- out.lda.lst[[jmark]]
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  topics.mat <- tm.result$topics
  dat.umap <- DoUmapAndLouvain(topics.mat, jsettings)
  if (jmark != "H3K4me1-H3K9me3"){
    dat.metas.tmp <- subset(dat.metas[[jmark]], select = c(-umap1, -umap2, -louvain))
    dat.umap.annot <- dat.umap %>%
      left_join(., dat.metas.tmp) %>%
      mutate(louvain = cluster)
  } else {
    dat.umap.annot <- dat.umap  %>%
      mutate(cluster = louvain)
  }
  return(dat.umap.annot)
})


# Plot UMAP  --------------------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.lst <- lapply(jmarks.withdbl, function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)


# Plot K4e1 and K9me3 proper colors ---------------------------------------

m.lst <- lapply(c("H3K4me1", "H3K9me3"), function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point() + 
    # scale_color_identity() + 
    scale_color_identity(guide = "legend",
                         labels = unique(dat.umap.lst[[jmark]]$cluster), 
                         breaks = unique(dat.umap.lst[[jmark]]$clustercol)) + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], cols = 2)

m.lst <- lapply(c("H3K4me1", "H3K9me3"), function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = stypecol)) + 
    geom_point() + 
    # scale_color_identity() + 
    scale_color_identity(guide = "legend",
                         labels = unique(dat.umap.lst[[jmark]]$batch), 
                         breaks = unique(dat.umap.lst[[jmark]]$stypecol)) + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], cols = 2)



# 
# 
# # Write outlda and cluster  ---------------------------------------------
# 
# for (jmark in jmarks.withdbl){
#   print(jmark)
#   jmarkout <- ifelse(jmark == "H3K4me1-H3K9me3", "H3K4me1xH3K9me3", jmark)
#   fnameout <- paste0("ClusterAnnot.lda_and_datmerged.k4_k9_", jsuffix2, ".", jmarkout, ".", Sys.Date(), ".RData")
#   outf <- file.path(outdir, fnameout)
#   out.lda <- out.lda.lst[[jmark]]
#   count.mat <- count.mat.lst[[jmark]]
#   dat.merge <- dat.umap.lst[[jmark]]
#   print(dim(count.mat))
#   save(out.lda, count.mat, dat.merge, file = outf)
# }
# 
# lapply(out.lda.lst, function(jout) length(jout@terms))
# lapply(out.lda.lst, function(jout) length(jout@documents))
# 
# 
# 
