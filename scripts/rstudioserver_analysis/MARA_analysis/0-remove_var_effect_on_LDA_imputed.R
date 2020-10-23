# Jake Yeung
# Date of Creation: 2020-08-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/0-remove_var_effect_on_LDA_imputed.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks



# Load LDA inputs ---------------------------------------------------------

jmark.ref <- "H3K4me1"
ntopics <- 30
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load bins ---------------------------------------------------------------

inf.bin <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000.imputevarfilt.lessstringent/lda_outputs.counts_table_var_filt.", jmark.ref, ".imputevar_0.75.K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark.ref, ".imputevar_0.75.K-30.Robj")
load(inf.bin, v=T)

out.lda.bins <- out.lda
tm.result.bins <- posterior(out.lda.bins)

dat.umap.bins <- DoUmapAndLouvain(tm.result.bins$topics, jsettings)

jchromos <- paste("chr", seq(25), sep = "")

dat.impute.log <- t(log2(tm.result.bins$topics %*% tm.result.bins$terms))
dat.var <- CalculateVarAll(dat.impute.log, jchromos)


# Load peaks --------------------------------------------------------------




hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks/lda_outputs.", jmark.ref, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-", ntopics, ".binarize.FALSE/ldaOut.", jmark.ref, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-", ntopics, ".Robj"))
load(inf, v=T)

# load inf bins for calculating variance


tm.result <- posterior(out.lda)


dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)


# Load annots -------------------------------------------------------------

dat.annot.lst <- lapply(jmarks, function(jmark){
  inf.annot <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt"))
  assertthat::assert_that(file.exists(inf.annot))
  dat.annot <- fread(inf.annot)
  dat.annot$mark <- jmark
  return(dat.annot)
})

dat.annot.mark <- subset(dat.annot.lst$H3K4me1, select = c(cell, cluster, plate)) %>%
  left_join(., dat.var)

dat.merge <- left_join(dat.umap, dat.annot.mark)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point () +
  theme_bw() +  
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(jmark.ref, ntopics))


ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point () +
  theme_bw() +  
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(jmark.ref, ntopics))


# does PCA solve this? 


# load exprs mat ----------------------------------------------------------

inf.exprs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/H3K4me1/mara_input/count_mats_peaks_norm/ldaOut.H3K4me1.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.keepNbins_250.txt"

exprs <- read.table.handlerows(inf.exprs)

pca.out <- irlba(A = t(exprs))
# pca.out <- prcomp(t(exprs), center = TRUE, scale. = TRUE)

dat.pca <- data.frame(cell = colnames(exprs), pca.out$u, stringsAsFactors = FALSE) %>%
  left_join(., dat.annot.mark %>% ungroup() %>% mutate(cell = make.names(cell)))

ggplot(dat.pca, aes(x = X1, y = X2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

