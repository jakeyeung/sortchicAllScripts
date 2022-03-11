# Jake Yeung
# Date of Creation: 2022-01-07
# File: ~/projects/scchic/scripts/macbook_analysis_2021/integrated_analysis_again/replot_integrated_data.R
# Replot on macbook because Rstudio server died

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

# Load data  --------------------------------------------------------------

inf.baccin <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/analysis_scrnaseq_atacseq/Baccin_scRNAseq_bonemarrow_no_niche.2021-09-05.RData"
load(inf.baccin, v=T)  # Blood, dat.umap
dat.meta.baccin <- dat.umap

inf.cusanovich.lda <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/analysis_scrnaseq_atacseq/scATACseq_Cusanovich_10kbTSS_LDA.rds"
dat.cusanovich.lda <- readRDS(inf.cusanovich.lda)

inf.meta.cusanovich <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/analysis_scrnaseq_atacseq/scATACseq_Cusanovich_meta_withumap.txt"
dat.meta.cusanovich <- fread(inf.meta.cusanovich)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597", "#eb9d01", "#7fbedf")
cbPalette2 <- c(cbPalette, "#f09de6")

ggplot(dat.meta.baccin, aes(x = UMAP_1, y = UMAP_2, color = celltype)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Integrated analysis  ----------------------------------------------------

inf.baccin.cca <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/analysis_scrnaseq_atacseq/scRNAseqBaccin_sortChICZeller_CCA_comparison.RData"
load(inf.baccin.cca, v=T)

inf.cusanovich.cca <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/analysis_scrnaseq_atacseq/scATACseqCusanovich_sortChICZeller_CCA_comparison.RData"
load(inf.cusanovich.cca, v=T)
dat.umap.atac.chic.annot <- dat.umap.cca.annot



# plot outputs ------------------------------------------------------------

head(dat.umap.cca.rna.chic.annot)
experis <- unique(dat.umap.cca.rna.chic.annot$experi); names(experis) <- experis

m.lst <- lapply(experis, function(jexperi){
  m <- ggplot(dat.umap.cca.rna.chic.annot %>% filter(experi == jexperi), aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jexperi) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})


# label individual clusters
jclst <- "small pre-B."
jclst <- "pro-B"
jclst <- "large pre-B."
jclst <- "Ery/Mk prog."

jclsts <- sort(unique(dat.umap.cca.rna.chic.annot$cluster))

outdir <- "/Users/yeung/data/scchic/from_cluster_2021/analysis_scrnaseq_atacseq/plots_from_macbookpro"
outpdf.rna <- file.path(outdir, paste0("check_cluster_cca_rna.", Sys.Date(), ".pdf"))
outpdf.atac <- file.path(outdir, paste0("check_cluster_cca_atac.", Sys.Date(), ".pdf"))


jexperi <- "RNA"
jclsts.rna <- sort(unique(subset(dat.umap.cca.rna.chic.annot, experi == jexperi)$cluster))

# pdf(outpdf.rna, useDingbats = FALSE)
JFuncs::multiplot(m.lst[[1]], m.lst[[2]])
print(m.lst[[1]])
print(m.lst[[2]])
# for (jclst in jclsts.rna){
#   print(jclst)
#   m.check <- ggplot(dat.umap.cca.rna.chic.annot %>% filter(experi == jexperi) %>% mutate(cluster = ifelse(cluster == jclst, jclst, "zNot")) %>% arrange(desc(cluster)), 
#                     aes(x = umap1, y = umap2, color = cluster)) + 
#     geom_point() + 
#     scale_color_manual(values = cbPalette) + 
#     ggtitle(paste(jexperi, jclst)) + 
#     theme_bw() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   print(m.check)
# }
# dev.off()

# merge certain cluster together

# Do for ATAC -------------------------------------------------------------

# head(dat.umap.cca.rna.chic.annot)
experis <- unique(dat.umap.atac.chic.annot$experi); names(experis) <- experis

print(experis)
m.lst <- lapply(experis, function(jexperi){
  m <- ggplot(dat.umap.atac.chic.annot %>% filter(experi == jexperi), aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette2) + 
    ggtitle(jexperi) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})


# label individual clusters
jclst <- "Activated B cells"
jclst <- "B cells"
jclst <- "Dendritic cells"
jclst <- "NK cells"
jclst <- "T cells"
jclst <- "Hematopoietic progenitors"


jexperi <- "ATAC"
# jclsts.atac <- sort(unique(subset(dat.umap.atac.chic.annot, experi == jexperi)$cluster))
# pdf(outpdf.atac, useDingbats = FALSE)
# JFuncs::multiplot(m.lst[[1]], m.lst[[2]])
# print(m.lst[[1]])
# print(m.lst[[2]])
# for (jclst in jclsts.atac){
#   print(jclst)
#   m.check <- ggplot(dat.umap.atac.chic.annot %>% filter(experi == jexperi) %>% 
#                       mutate(cluster = ifelse(cluster == jclst, jclst, "zNot")) %>% 
#                                arrange(desc(cluster)), 
#                     aes(x = umap1, y = umap2, color = cluster)) + 
#     geom_point() + 
#     scale_color_manual(values = cbPalette) + 
#     ggtitle(paste(jexperi, jclst)) + 
#     theme_bw() + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   print(m.check)
# }
# dev.off()



