# Jake Yeung
# Date of Creation: 2020-02-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_from_hiddendomains/2-downstream_GLMPCA.R
# GLMPCA

rm(list=ls())

jstart <- Sys.time()

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
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load LDA and GLMPCA -----------------------------------------------------

# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
ntopics <- "30"
jbinskeep <- 250

for (jmark in jmarks){
  
  
  print(jmark)
  
  inf.glmpca <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs_fromhiddendomains.NoVarCorrection.KeepBestPlates2/PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.mergebinsize_1000.binskeep_", jbinskeep, ".covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.ntopics_", ntopics, ".2020-02-11.RData")
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains/lda_outputs.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".binarize.FALSE/ldaOut.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".Robj")
  inf.lda.bins <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
  
  assertthat::assert_that(file.exists(inf.glmpca))
  assertthat::assert_that(file.exists(inf.lda))
  assertthat::assert_that(file.exists(inf.lda.bins))
  
  assertthat::assert_that(endsWith(inf.glmpca, suffix = ".RData"))
  jbase <- paste0("PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_", jbinskeep, ".ntopics_", ntopics, ".2020-02-11")
  
  # outname <- paste0("countmat_", gsub(".RData$", ".txt", basename(inf.glmpca)))
  outname <- paste0("countmat_", jbase, ".txt")
  outnamepdf <- paste0("countmat_", jbase, ".pdf")
  outdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_input/count_mats_peaks_norm")
  dir.create(outdir, recursive = TRUE)
  outf <- file.path(outdir, outname)
  outpdf <- file.path(outdir, outnamepdf)
  
  # if (file.exists(outf)){
  #   print(paste("Outf exists,",  outf, " skipping"))
  #   next
  # }
  
  
  load(inf.glmpca, v=T)
  
  load(inf.lda, v=T)
  out.lda.hd <- out.lda
  count.mat.hd <- count.mat
  
  load(inf.lda.bins, v=T)
  out.lda.bins <- out.lda
  count.mat.bins <- count.mat
  
  pdf(outpdf, useDingbats = FALSE)
  
  dat.umap.hd.glmpca <- DoUmapAndLouvain(glm.out$factors, jsettings) %>%
    dplyr::rename(louvain.hd.glmpca = louvain) %>% 
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"),
           cond = GetCondFromSamp(as.character(cell), mark = jmark))
  
  dat.umap.hd.lda <- DoUmapAndLouvain(posterior(out.lda.hd)$topics, jsettings) %>%
    dplyr::rename(louvain.hd.lda = louvain) %>% 
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"),
           cond = GetCondFromSamp(as.character(cell), mark = jmark))
  
  dat.umap.bins.lda <- DoUmapAndLouvain(posterior(out.lda.bins)$topics, jsettings) %>%
    dplyr::rename(louvain.bins.lda = louvain) %>% 
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"),
           cond = GetCondFromSamp(as.character(cell), mark = jmark))
  
  ggplot(dat.umap.hd.glmpca, aes(x = umap1, y = umap2, color = louvain.hd.glmpca)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle("HiddenDomains GLMPCA")
  ggplot(dat.umap.hd.glmpca, aes(x = umap1, y = umap2, color = louvain.hd.glmpca)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle("HiddenDomains GLMPCA") + facet_wrap(~cond)
  
  ggplot(dat.umap.hd.lda, aes(x = umap1, y = umap2, color = louvain.hd.lda)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle("HiddenDomains LDA")
  ggplot(dat.umap.hd.lda, aes(x = umap1, y = umap2, color = louvain.hd.lda)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle("HiddenDomains LDA") + facet_wrap(~cond)
  
  ggplot(dat.umap.bins.lda, aes(x = umap1, y = umap2, color = louvain.bins.lda)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle("Bins LDA")
  ggplot(dat.umap.bins.lda, aes(x = umap1, y = umap2, color = louvain.bins.lda)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
    ggtitle("Bins LDA") + facet_wrap(~cond)
  
  
  # plot variance
  
  dat.umap.bins.lda.var <- left_join(dat.umap.bins.lda, subset(dat.merge2, select = c(cell, ncuts, ncuts.var.log2.CenteredAndScaled)))
  dat.umap.hd.lda.var <- left_join(dat.umap.hd.lda, subset(dat.merge2, select = c(cell, ncuts, ncuts.var.log2.CenteredAndScaled)))
  dat.umap.hd.glmpca.var <- left_join(dat.umap.hd.glmpca, subset(dat.merge2, select = c(cell, ncuts, ncuts.var.log2.CenteredAndScaled)))
  
  ggplot(dat.umap.bins.lda.var, aes(x = umap1, y = umap2, color = ncuts.var.log2.CenteredAndScaled)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1)  + ggtitle("LDA bins, var from bins")
  ggplot(dat.umap.hd.lda.var, aes(x = umap1, y = umap2, color = ncuts.var.log2.CenteredAndScaled)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + 
    facet_wrap(~cond) + ggtitle('LDA hd, var from bins')
  
  ggplot(dat.umap.hd.glmpca.var, aes(x = umap1, y = umap2, color = ncuts.var.log2.CenteredAndScaled)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + 
    ggtitle('GLMPCA hd, var from bins')
  ggplot(dat.umap.hd.glmpca.var, aes(x = umap1, y = umap2, color = ncuts.var.log2.CenteredAndScaled)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + 
    facet_wrap(~cond) + ggtitle('GLMPCA hd, var from bins')
  
  # Recalculalte variance ---------------------------------------------------
  
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  tm.result <- posterior(out.lda.hd)
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.var.hd <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.umap.hd.lda.var2 <- left_join(dat.umap.hd.lda, dat.var.hd)
  dat.umap.hd.glmpca.var2 <- left_join(dat.umap.hd.glmpca, dat.var.hd)
  
  ggplot(dat.umap.hd.lda.var2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle("HiddenDomains LDA, var from peaks")
  ggplot(dat.umap.hd.lda.var2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle("HiddenDomains LDA, var from peaks") + 
    facet_wrap(~cond)
  
  ggplot(dat.umap.hd.glmpca.var2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + ggtitle("HiddenDomains GLMPCA, var from peaks") 
  ggplot(dat.umap.hd.glmpca.var2, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) +  ggtitle("HiddenDomains GLMPCA, var from peaks")
    facet_wrap(~plate)
  
  # check if PCA associated with variance?
  
  pca.out <- prcomp(tm.result$topics, center = TRUE, scale. = TRUE)
  dat.pca <- data.frame(cell = rownames(pca.out$x), pc1 = pca.out$x[, 1], pc2 = pca.out$x[, 2], stringsAsFactors = FALSE) %>%
    left_join(., dat.var.hd) %>%
    left_join(., subset(dat.merge2, select = c(cell, ncuts.var.log2.CenteredAndScaled)))
  
  ggplot(dat.pca, aes(x = pc1, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point()  + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_viridis_c(direction = -1) + 
    ggtitle("HiddenDomains PCA var recalculated from hiddendomains")
  
  ggplot(dat.pca, aes(x = pc1, y = ncuts.var.log2.CenteredAndScaled, color = cell.var.within.sum.norm)) + geom_point()  + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_viridis_c(direction = -1) + 
    ggtitle("HiddenDomains PCA var recalculated from hiddendomains")
  
  ggplot(dat.pca, aes(x = pc2, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point()  + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_viridis_c(direction = -1) + 
    ggtitle("HiddenDomains PCA var recalculated from hiddendomains")
  
  ggplot(dat.pca, aes(x = pc2, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point()  + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_viridis_c(direction = -1) + 
    ggtitle("HiddenDomains PCA var recalculated from hiddendomains")
  
  
  # Write output to file for MARA -------------------------------------------
  
  tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = as.matrix(t(glm.out$loadings)))
  dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)
  rnames <- rownames(dat.impute.glm)  # chr:start-end;gene
  # strip off gene
  rnames.strip <- sapply(rnames, function(x) strsplit(x, ";")[[1]][[1]])
  dat.impute.glm.output <- data.frame(Gene.ID = rnames.strip, dat.impute.glm, stringsAsFactors = FALSE)
  
  
  fwrite(dat.impute.glm.output, file = outf, row.names = FALSE, sep = "\t")
  
  print("Done")
  print(dat.impute.glm.output[1:5, 1:5])
  
  print(Sys.time() - jstart)
  dev.off()
  
}


