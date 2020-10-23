# Jake Yeung
# Date of Creation: 2020-09-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/LDA_LSI_glmpca_downstream.R
# Downstream of LDA, LSI, and GLMPCA 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(irlba)

library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Get Hoesch information --------------------------------------------------



# Load spikeins -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.rds"
dat.spikeins <- readRDS(inf.spikeins)

dat.spikeins <- AddCellCycleLabel.bydat(dat.spikeins)

# Load LDA  ---------------------------------------------------------------

jmark <- jmarks[[3]]
# jsuffix <- "cellcyclefilt"
jsuffix <- "AllMerged"
jsuffix <- "G1filt"
jtopn <- 5000

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

if (jsuffix == "AllMerged"){
  allmerged <- TRUE
} else {
  allmerged <- FALSE
}

if (!allmerged){
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562.round2"
} else {
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562.round2.AllMerged"
}
dir.create(outdir)


for (jmark in jmarks){
  jtitle <- paste(jmark, jsuffix, jtopn)
  print(jtitle)
  
  outpdf <- file.path(outdir, paste0("LDA_GLMPCA_LSI_downstream.", jmark, ".", jsuffix, ".", jtopn, ".pdf"))
  
  pdf(outpdf, useDingbats = FALSE)
  
  # for cellcycle?
  if (!allmerged){
    inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_K562_spikein_round2.genes_filt/lda_outputs.K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".K-30.binarize.FALSE/ldaOut.K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".K-30.Robj")
  } else {
    inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_K562_spikein_round2.merged.genes_filt/lda_outputs.K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".K-30.binarize.FALSE/ldaOut.K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".K-30.Robj")
  }
  # for all merged
  assertthat::assert_that(file.exists(inf.lda))
  
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  
  dat.var.raw <- CalculateVarRaw(count.mat)
  
  jchromos <- paste("chr", seq(19), sep = "")
  
  dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))
  
  dat.var.imputed <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.var <- left_join(dat.var.raw, dat.var.imputed, by = "cell")
  
  # Load glmpca -------------------------------------------------------------
  
  # inf.glmpca <- paste0("/home/jyeu  ng/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2/K562_count_tables_50000.", jmark, ".", jsuffix, ".glmpcaout.penalty_1.maxiter_10000.RData")
  if (!allmerged){
    inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2.genes_filt/K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".glmpcaout.penalty_1.maxiter_10000.RData")
  } else {
    inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2.merged.genes_filt/K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".glmpcaout.penalty_1.maxiter_10000.RData")
  }
  assertthat::assert_that(file.exists(inf.glmpca))
  load(inf.glmpca, v=T)
  
  
  # Plot outputs ------------------------------------------------------------
  
  
  dat.umap.lda <- DoUmapAndLouvain(tm.result$topics, jsettings)
  
  dat.umap.lda.annot <- left_join(dat.umap.lda, dat.var, by = "cell") %>%
    left_join(., dat.spikeins, by = c("cell" = "samp"))
  
  
  m <- ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.lda.annot, aes(x = umap1, y = umap2, color = log10(chromocounts))) + 
    geom_point() + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.lda.annot, aes(x = cell.var.within.sum.norm, y = ncuts.var)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.lda.annot, aes(x = cell.var.within.sum.norm, y = log2(ncuts.var))) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  
  
  # Plot GLMPCA -------------------------------------------------------------
  
  dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings)
  dat.umap.glmpca$dim1 <- glmpcaout$factors$dim1
  dat.umap.glmpca$dim2 <- glmpcaout$factors$dim2
  
  dat.umap.glmpca.annot <- dat.umap.glmpca %>%
    left_join(., dat.spikeins, by = c("cell" = "samp"))  %>%
    left_join(., dat.var, by = "cell")
  
  m <- ggplot(dat.umap.glmpca.annot, aes(x = dim1, y = dim2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() +  
    scale_color_viridis_c(direction = 1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.glmpca.annot, aes(x = dim1, y = dim2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() +  
    scale_color_viridis_c(direction = 1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.glmpca.annot, aes(x = dim1, y = dim2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() +  
    scale_color_viridis_c(direction = 1) +  
    facet_wrap(~cellcycle.str) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.glmpca.annot, aes(x = dim1, y = dim2, color = log2(chromocounts))) + 
    geom_point() + 
    theme_bw() +  
    scale_color_viridis_c(direction = 1) +  
    facet_wrap(~cellcycle.str) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  m <- ggplot(dat.umap.glmpca.annot, aes(x = dim1, y = dim2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() +  
    scale_color_viridis_c(direction = 1) +  
    facet_wrap(~cellcycle.str) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m)
  
  
  
  # Do LSI ------------------------------------------------------------------
  
  lsi.out <- RunLSI(as.matrix(count.mat))
  
  dat.lsi.umap <- DoUmapAndLouvain(lsi.out$u, jsettings) %>%
    left_join(., dat.var, by = "cell") %>%
    left_join(., dat.spikeins, by = c("cell" = "samp"))
  
  m <- ggplot(dat.lsi.umap, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.lsi.umap, aes(x = umap1, y = umap2, color = log2(chromocounts))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.lsi.umap, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  dev.off()
  
}


