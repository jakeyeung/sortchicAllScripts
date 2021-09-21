# Jake Yeung
# Date of Creation: 2021-08-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/14-compare_atacseq_vs_sortchic_genesets.R
# 

rm(list=ls())

library(scchicFuncs)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 7

hubprefix <- "/home/jyeung/hub_oudenaarden"

inmain <- "jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq"

# Load atacseq -------------------------------------------------------------


inf.tm.atac <- file.path(hubprefix, inmain, "scATACseq_Cusanovich_10kbTSS_LDA.rds")
inf.mat.atac <- file.path(hubprefix, inmain, "scATACseq_Cusanovich_10kbTSS_countmat.rds")
inf.meta.atac <- file.path(hubprefix, inmain, "scATACseq_Cusanovich_meta.txt")

tm.atac <- readRDS(inf.tm.atac)
mat.atac <- readRDS(inf.mat.atac)
meta.atac <- fread(inf.meta.atac)


# Load sortchic -----------------------------------------------------------

inf.tm.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_10kbTSS_LDA.rds")
inf.mat.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_10kbTSS_countmat.rds")
inf.meta.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_meta.txt")

tm.chic <- readRDS(inf.tm.chic)
mat.chic <- readRDS(inf.mat.chic)
meta.chic <- fread(inf.meta.chic)


# Make UMAPs --------------------------------------------------------------

dat.umap.atac <- DoUmapAndLouvain(tm.atac$topics, jsettings = jsettings)
dat.umap.chic <- DoUmapAndLouvain(tm.chic$topics, jsettings = jsettings)

dat.umap.atac.annot <- left_join(dat.umap.atac, meta.atac)
dat.umap.chic.annot <- left_join(dat.umap.chic, meta.chic)


# impute
dat.imputed.atac <- t(log2(tm.atac$topics %*% tm.atac$terms))
dat.imputed.chic <- t(log2(tm.chic$topics %*% tm.chic$terms))

m.atac <- ggplot(dat.umap.atac.annot, aes(x = umap1, y = umap2, color = cell_label)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  ggtitle("scATACseq Cusanovich +/-5 TSS bins") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.chic <- ggplot(dat.umap.chic.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("sortChIC Zeller +/-5 TSS bins") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load gene sets ----------------------------------------------------------

inf.gsets <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/post_submission/genesets/gene_sets_from_sortChIC.txt")
dat.gsets <- fread(inf.gsets)



# Show gene sets onto UMAP  -----------------------------------------------

# do eryths
jjset.vec <- unique(dat.gsets$jset)

# jjset <- "Eryths"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/plots/scATACseq_sortChIC_comparison.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)
print(m.atac)
print(m.atac + facet_wrap(~cell_label))
print(m.chic)
print(m.chic + facet_wrap(~cluster))
for (jjset in jjset.vec){
  
  jregions <- subset(dat.gsets, jset == jjset)$gene
  
  dat.imputed.atac.filt <- dat.imputed.atac[jregions, ]
  assertthat::assert_that(nrow(dat.imputed.atac.filt) > 0)
  exprs.atac <- colMeans(dat.imputed.atac.filt)
  dat.exprs.atac <- data.frame(cell = names(exprs.atac), gexprs = exprs.atac, stringsAsFactors = FALSE) %>%
    left_join(., dat.umap.atac.annot)
  
  dat.imputed.chic.filt <- dat.imputed.chic[jregions, ]
  assertthat::assert_that(nrow(dat.imputed.chic.filt) > 0)
  exprs.chic <- colMeans(dat.imputed.chic.filt)
  dat.exprs.chic <- data.frame(cell = names(exprs.chic), gexprs = exprs.chic, stringsAsFactors = FALSE) %>%
    left_join(., dat.umap.chic.annot)
  
  
  jtitle <- paste("average levels +/-5 kb around", jjset, "genes")
  
  m.exprs.atac <- ggplot(dat.exprs.atac, aes(x = umap1, y = umap2, color = gexprs)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtitle, "ATAC") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c() 
  
  m.exprs.chic <- ggplot(dat.exprs.chic, aes(x = umap1, y = umap2, color = gexprs)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtitle, "ChIC") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c() 
  
  JFuncs::multiplot(m.exprs.atac, m.exprs.chic, cols = 2)
  
}


# Check S100a8
# jgene <- "S100a8"
jgenes <- c("S100a8", "Cebpe", "Ptpre", "S100a2", "S100a9")
for (jgene in jgenes){
  
  # jgene <- "S100a8"
  jregions <- subset(dat.gsets, symbol == jgene)$gene
  
  dat.imputed.atac.filt <- dat.imputed.atac[jregions, ]
  if (length(jregions) > 1){
    exprs.atac <- colMeans(dat.imputed.atac.filt)
  } else {
    exprs.atac <- dat.imputed.atac.filt
  }
  dat.exprs.atac <- data.frame(cell = names(exprs.atac), gexprs = exprs.atac, stringsAsFactors = FALSE) %>%
    left_join(., dat.umap.atac.annot)
  
  
  
  
  dat.imputed.chic.filt <- dat.imputed.chic[jregions, ]
  if (length(jregions) > 1){
    exprs.chic <- colMeans(dat.imputed.chic.filt)
  } else {
    exprs.chic <- dat.imputed.chic.filt
  }
  # exprs.chic <- colMeans(dat.imputed.chic.filt)
  dat.exprs.chic <- data.frame(cell = names(exprs.chic), gexprs = exprs.chic, stringsAsFactors = FALSE) %>%
    left_join(., dat.umap.chic.annot)
  
  
  jtitle <- paste("average levels +/-5 kb around", jgene)
  
  m.exprs.atac <- ggplot(dat.exprs.atac, aes(x = umap1, y = umap2, color = gexprs)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtitle, "ATAC") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c() 
  
  m.exprs.chic <- ggplot(dat.exprs.chic, aes(x = umap1, y = umap2, color = gexprs)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtitle, "ChIC") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c() 
  
  JFuncs::multiplot(m.exprs.atac, m.exprs.chic, cols = 2)
  
}
dev.off()

# # write UMAP to meta
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq"
# outf.meta.atac <- file.path(outdir, paste0("scATACseq_Cusanovich_meta_withumap.txt"))
# outf.meta.chic <- file.path(outdir, paste0("H3K4me1-sortChIC_Zeller_meta_withumap.txt"))
# 
# fwrite(dat.umap.atac.annot, file = outf.meta.atac, sep = "\t")
# fwrite(dat.umap.chic.annot, file = outf.meta.chic, sep = "\t")
# 
