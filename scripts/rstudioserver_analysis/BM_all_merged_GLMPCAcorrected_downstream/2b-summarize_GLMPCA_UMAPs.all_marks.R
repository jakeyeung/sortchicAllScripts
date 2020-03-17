# Jake Yeung
# Date of Creation: 2020-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/2-summarize_GLMPCA_UMAPs.R
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


# Constants ---------------------------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

# set conditions ----------------------------------------------------------

mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jconds <- c("Unenriched", "AllMerged")

# jmark <- "H3K4me1"
# jexperi <- "Unenriched"

jcondsmarks <- levels(interaction(jconds, jmarks, sep = "_"))

pdfdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMLPCA_outputs.KeepBestPlates2.downstream.showBeforeAfter"
outname <- paste0("GLMPCA_downstream_alllconds.mergesize_", mergesize, 
                  ".nbins_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".pdf")
# pdf(file.path(pdfdir, paste0(basename(inf.glm), ".pdf")))
pdf(file.path(pdfdir, outname), useDingbats = FALSE)

for (jcondmark in jcondsmarks){
  print(jcondmark)
  jexperi <- strsplit(jcondmark, "_")[[1]][[1]] 
  jmark <- strsplit(jcondmark, "_")[[1]][[2]]
  
  inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
  assertthat::assert_that(file.exists(inf.glm))
  
  
  # Load LDA output ---------------------------------------------------------
  
  
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.Robj")
  assertthat::assert_that(file.exists(inf.lda))
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  # label  tpic names
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "")
  
  
  dat.umap.lda <- DoUmapAndLouvain(topics.mat = posterior(out.lda)$topics, jsettings = jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"))
  
  # Load GLMPCA correctoin --------------------------------------------------
  
  
  load(inf.glm, v=T)
  
  tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = t(as.matrix(glm.out$loadings)))
  dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)
  
  # Pload umap before and after  ---------------------------------------------
  
  
  
  dat.umap.glm <- DoUmapAndLouvain(tm.result.glm$topics, jsettings = jsettings)
  dat.umap.glm$louvain.glm <- dat.umap.glm$louvain; dat.umap.glm$louvain <- NULL
  dat.umap.glm <- left_join(dat.umap.glm, subset(dat.umap.lda, select = c(cell, louvain, plate)))
  
  
  m.lda <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(jcondmark, "From LDA")
  m.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(jcondmark, "GLMPCA, colored by LDA louvain")
  m.glm.new.louvain <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(jcondmark, "GLMPCA new louvain")
  
  # plot different plates
  m.lda.plates <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(jcondmark, "From LDA") + 
    facet_wrap(~plate)
  m.glm.plates <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(jcondmark, "GLMPCA, colored by LDA louvain") + 
    facet_wrap(~plate)
  m.glm.new.louvain.plates <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + ggtitle(jcondmark, "GLMPCA new louvain") + 
    facet_wrap(~plate)
  
  multiplot(m.lda, m.glm, m.glm.new.louvain, cols = 3)
  lapply(list(m.lda.plates, m.glm.plates, m.glm.new.louvain.plates), print)
}
dev.off()

