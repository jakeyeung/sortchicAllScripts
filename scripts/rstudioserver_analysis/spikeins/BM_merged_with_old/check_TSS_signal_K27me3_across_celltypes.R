# Jake Yeung
# Date of Creation: 2020-11-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/check_TSS_signal_K27me3_across_celltypes.R
# 


rm(list=ls())

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(hash)
library(igraph)
library(umap)
library(ggrepel)

library(scchicFuncs)

jstart <- Sys.time()

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load LDA (contains countmat)  ---------------------------------------------------------------

# ncores <- 8
hubprefix <- "/home/jyeung/hub_oudenaarden"
# jtype <- "hiddendomains"
jtype <- "TSS"
# jdist <- "TES"
jdist <- 10000


# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K27me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K9me3"); names(jmarks) <- jmarks
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/global_hist_mod_BM"
for (jmark in jmarks){
  print(jmark)
  
  outpdf <- file.path(outdir, paste0("BM_signal_distributions_genomewide.", jtype, ".dist_", jdist, ".", jmark, ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  inf.spikein <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.blfix/spikein_info_BM_round2_all.blfix.txt")
  
  dat.spikein.all <- fread(inf.spikein) %>%
    mutate(cell = samp)
  
  indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04")
  inf.annot <- file.path(indir.annot, paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt"))
  dat.annot <- fread(inf.annot)
  
  indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_", jtype, ".dist_", jdist))
  fname <- paste0("lda_outputs.count_mat_from_", jtype, ".", jmark, ".dist_", jdist, ".K-30.binarize.FALSE/ldaOut.count_mat_from_", jtype, ".", jmark, ".dist_", jdist, ".K-30.Robj")
  
  inf.tmp <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.tmp))
  
  load(inf.tmp, v=T)
  
  dat.spikein <- subset(dat.spikein.all, cell %in% colnames(count.mat))
  
  dat.merge <- left_join(dat.spikein, dat.annot %>% dplyr::select(c('cell', 'cluster', 'plate')))
  
  cnames.keep.lst <- split(x = dat.annot$cell, f = dat.annot$cluster)
  
  pseudobulk <- SumAcrossClusters(count.mat, cnames.keep.lst)
  
  m0 <- ggplot(dat.merge, aes(x = cluster, y = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    facet_wrap(~plate) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m0)
  
  ncells.per.clstr <- dat.merge %>%
    group_by(cluster) %>%
    summarise(ncells = length(cell))
  
  spikeins.per.clstr <- dat.merge %>%
    group_by(cluster) %>%
    summarise(spikeincounts = sum(spikeincounts))
  
  jmerge <- left_join(ncells.per.clstr, spikeins.per.clstr)
  
  m.ncells <- ggplot(jmerge, aes(x = ncells, y = spikeincounts, label = cluster)) + 
    geom_text_repel() + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    scale_x_log10() + 
    scale_y_log10() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.ncells)
  
  pseudobulk.mat <- do.call(cbind, pseudobulk)
  
  # assertthat::assert_that(identical(ncells.per.clstr$cluster, colnames(pseudobulk.mat)))
  # assertthat::assert_that(identical(spikeins.per.clstr$cluster, colnames(pseudobulk.mat)))
  # 
  # # pseudobulk.mat.cellnorm <- sweep(pseudobulk.mat, MARGIN = 2, STATS = ncells.per.clstr$ncells, FUN = "/")
  pseudobulk.mat.cellnorm <- sweep(pseudobulk.mat, MARGIN = 2, STATS = spikeins.per.clstr$spikeincounts, FUN = "/")
  
  pseudobulk.mat.cellnorm.long <- pseudobulk.mat.cellnorm %>%
    melt()
  colnames(pseudobulk.mat.cellnorm.long) <- c("jname", "ctype", "ratio")
  
  m1 <- ggplot(pseudobulk.mat.cellnorm.long, aes(x = log2(ratio), fill = ctype)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(~ctype, ncol = 1) + 
    ggtitle(jmark) + 
    theme_bw() + 
    scale_fill_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m2 <- ggplot(pseudobulk.mat.cellnorm.long, aes(x = log2(ratio), fill = ctype)) + 
    geom_density(alpha = 0.5) + 
    facet_wrap(~ctype) + 
    ggtitle(jmark) + 
    theme_bw() + 
    scale_fill_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  print(m1)
  print(m2)
  
  dev.off()
  
  
}

