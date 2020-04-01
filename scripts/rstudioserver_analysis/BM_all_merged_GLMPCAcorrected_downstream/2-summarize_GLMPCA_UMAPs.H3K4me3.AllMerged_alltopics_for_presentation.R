# Jake Yeung
# Date of Creation: 2020-02-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/2-summarize_GLMPCA_UMAPs.H3K4me1.AllMerged_alltopics.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)
library(igraph)
library(umap)

library(mixtools)

library(RColorBrewer)

# Constants ---------------------------------------------------------------

write.plots <- TRUE

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#ff9f7d", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#eb9d01", "#7fbedf", "#009E73")



jthres.frac.cells <- 0.9  # threshold for imputing NAs to celltype 
# (e.g. NA associarted to louvain cluster that are 90% of one celltype)

# topic14 is not convergent :( 
# jtopics <- paste("topic", c(11, 18, 10, 4, 29, 8, 14, 23, 27, 5, 16, 7, 22), sep = "")
# names(jtopics) <- c("ILC-RoraPlus", "ILC-PrkchPlus", "Bcells-Cd83", "Dendritic", "Bcells-Cd47", "Basophils", "Unknown-Igf1r", "Neutrophils", "Eryth", "HSCs-Ephb2", "pDendritic", "HSCs-Hlf", "HSCs-Anxa2")

# jtopics <- paste("topic", c(11, 18, 10, 4, 29, 8, 23, 27, 5, 16, 7, 22), sep = "")
# names(jtopics) <- c("ILC-RoraPlus", "ILC-PrkchPlus", "Bcells-Cd83", "Dendritic", "Bcells-Cd47", "Basophils", "Neutrophils", "Eryth", "HSCs-Ephb2", "pDendritic", "HSCs-Hlf", "HSCs-Anxa2")

jtopics <- paste("topic", seq(30), sep = "")
names(jtopics) <- jtopics

if (length(jtopics) > length(cbPalette)){
  print('using auto colors...')
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  cbPalette <- col_vector[1:length(jtopics)]
}


names(jtopics) <- paste(names(jtopics), jtopics, sep = "_")
names.final <- names(jtopics)

# add jtopics at the end for easy reference later?
# names(jtopics) <- paste(names(jtopics), jtopics, sep = "_")
jnames.hash <- hash::hash(names(jtopics), names.final)  # allows clusters to be merged afterwards by having the same final name

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# jmark <- "H3K4me3"


pdfdir <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/primetime")
dir.create(pdfdir)
pdfname <- paste0("GLMPCA_vs_LDA.allmarks.", Sys.Date(), ".pdf")
pdfout <- file.path(pdfdir, pdfname)

if (write.plots){
  pdf(pdfout, width = 1440/72, height = 815/72, useDingbats = FALSE)
}
  
for (jmark in jmarks){
  jexperi <- "AllMerged"
  mergesize <- "1000"
  nbins <- "1000"
  jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
  jpenalty <- 1
  inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
  assertthat::assert_that(file.exists(inf.glm))
  
  jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
  assertthat::assert_that(file.exists(jinf.tss))
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  # Load LDA output ---------------------------------------------------------
  
  
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.Robj")
  assertthat::assert_that(file.exists(inf.lda))
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  # label  tpic names
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "")
  
  
  
  dat.umap.lda <- DoUmapAndLouvain(topics.mat = posterior(out.lda)$topics, jsettings = jsettings)
  
  # Load GLMPCA correctoin --------------------------------------------------
  
  
  load(inf.glm, v=T)
  
  sf <- data.frame(cell = rownames(glm.inits$X.mat), rawvar = glm.inits$X.mat[, 1], stringsAsFactors = FALSE)
 
  
  tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = t(as.matrix(glm.out$loadings)))
  dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)
  
  
  dat.umap.glm <- DoUmapAndLouvain(tm.result.glm$topics, jsettings = jsettings)
  dat.umap.glm$louvain.glm <- dat.umap.glm$louvain; dat.umap.glm$louvain <- NULL
  dat.umap.glm <- left_join(dat.umap.glm, subset(dat.umap.lda, select = c(cell, louvain)))
  
  # add variance
  dat.umap.lda <- left_join(dat.umap.lda, sf)
  dat.umap.glm <- left_join(dat.umap.glm, sf)
  
  
  m.lda <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, "From LDA"))
  m.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, "GLMPCA, colored by LDA louvain"))
  m.glm.new.louvain <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, "GLMPCA new louvain"))
  multiplot(m.lda, m.glm, m.glm.new.louvain, cols = 3)
  
  m.lda.var <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = rawvar)) + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, "From LDA"))
  m.glm.var <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = rawvar)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, "GLMPCA, colored by LDA louvain"))
  m.glm.new.louvain.var <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = rawvar)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, "GLMPCA new louvain"))
  multiplot(m.lda.var, m.glm.var, m.glm.new.louvain.var, cols = 3)
  
  # print one just for the legend
  m.legend <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = rawvar)) + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle(paste(jmark, "From LDA"))
  print(m.legend)
  
  
  
}


if (write.plots){
  dev.off()
}

