# Jake Yeung
# Date of Creation: 2021-01-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/playground/compare_LDA_with_PCA.R
# Quick check


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

library(glmpca)

library(irlba)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$spread <- 8
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")




# Load .metasmetas  -------------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})


# Load LDA output ---------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks.k9 <- "H3K9me3"  # add this separately

hubprefix <- "/home/jyeung/hub_oudenaarden"


# jmark <- "H3K9me3"
niter <- "100"
inf.glmpca.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/BM_from_matadj/glmpca.", jmark, ".from_matadj.platename_jrep.szname_none.niter_", niter, ".RData"))
    assertthat::assert_that(file.exists(inf.glmpca))
  } else {
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_plate.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData"))
  }
  return(inf.glmpca)
})


outs.lst <- lapply(inf.glmpca.lst, function(inf){
  load(inf, v=T)
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  return(list(dat.umap = dat.umap, count.mat = glm.inits$Y.filt))
})

# add metadata

dat.lda.lst <- lapply(jmarks, function(jmark){
  left_join(outs.lst[[jmark]]$dat.umap, subset(dat.metas[[jmark]], select = c(cell, cluster, clustercol)))
})


ggplot(dat.lda.lst$H3K9me3, aes(x = -1 * umap1, y = 1 * umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Do PCA ------------------------------------------------------------------

dat.lsi.lst <- lapply(jmarks, function(jmark){
  jout <- outs.lst[[jmark]]
  out.lsi <- RunLSI(jout$count.mat)
  dat.lsi <- DoUmapAndLouvain(out.lsi$u, jsettings = jsettings) %>%
    left_join(., subset(dat.metas[[jmark]], select = c(cell, cluster, clustercol)))
  return(dat.lsi)
})

dat.lsi.lst <- lapply(jmarks, function(jmark){
  jout <- outs.lst[[jmark]]
  out.lsi <- RunLSI(jout$count.mat)
  dat.lsi <- DoUmapAndLouvain(out.lsi$u, jsettings = jsettings) %>%
    left_join(., subset(dat.metas[[jmark]], select = c(cell, cluster, clustercol)))
  return(dat.lsi)
})

ggplot(dat.lsi.lst$H3K9me3, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.lsi.lst.test <- lapply(jmarks, function(jmark){
  jout <- outs.lst[[jmark]]
  out.lsi <- RunLSI(jout$count.mat, n.components = 30)
  dat.lsi <- DoUmapAndLouvain(out.lsi$u, jsettings = jsettings) %>%
    left_join(., subset(dat.metas[[jmark]], select = c(cell, cluster, clustercol)))
  return(dat.lsi)
})

# try standard PCA 

ggplot(dat.lsi.lst.test$H3K9me3, aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.metas$H3K9me3$cluster)) + 
    xlab("") + ylab("") + 
    theme_minimal() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

jmark <- "H3K9me3"
dat.svd.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jmat <- outs.lst[[jmark]]$count.mat
  jmat.norm <- log2(sweep(jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/") + 1 * 1000)
  jmat.norm2 <- scale(t(jmat.norm), center = TRUE, scale = TRUE)
  
  count.svd <- irlba(A = jmat.norm2, nv = 30, scale = FALSE, center = FALSE)
  
  rownames(count.svd$u) <- rownames(jmat.norm2)
  # jmat.norm2.eigengenes <- count.svd$u
  dat.svd <- DoUmapAndLouvain(count.svd$u, jsettings = jsettings)
  dat.svd.annot <- left_join(dat.svd, subset(dat.metas[[jmark]], select = c(cell, cluster, clustercol)))
  return(dat.svd.annot)
})

for (jmark in jmarks){
  m <- ggplot(dat.svd.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.metas$H3K9me3$cluster)) + 
    xlab("") + ylab("") + 
    theme_minimal() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}

# ggplot(dat.svd.annot, aes(x = umap1, y = umap2, color = clustercol)) + 
#   geom_point() + 
#   scale_color_manual(values = cbPalette) + 
#   scale_color_identity(guide = "legend", 
#                        labels = unique(dat.metas$H3K9me3$cluster)) + 
#   xlab("") + ylab("") + 
#   theme_minimal() + 
#   ggtitle(jmark) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 




# Compare LDA vs PCA ------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/for_presentations"
outpdf <- file.path(outdir, "lda_vs_pca.pdf")

pdf(file = outpdf, useDingbats = FALSE)
for(jmark in jmarks){
  m.lda <- ggplot(dat.lda.lst[[jmark]], aes(x = -1 * umap1, y = 1 * umap2, color = clustercol)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.metas$H3K9me3$cluster)) + 
    xlab("") + ylab("") + 
    theme_minimal() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  m.pca <- ggplot(dat.lsi.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.metas$H3K9me3$cluster)) + 
    theme_minimal() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  JFuncs::multiplot(m.lda, m.pca, cols = 2)
}
dev.off()





# Check other k9me3  ------------------------------------------------------


# inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows.cellfilt_binfilt/lda_outputs.count_mat.H3K9me3.match_dbl.cellfilt.binfilt.K-30.binarize.FALSE/ldaOut.count_mat.H3K9me3.match_dbl.cellfilt.binfilt.K-30.Robj")
inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_varfilt.blfix/lda_outputs.count_mat.H3K9me3.filt_0.15_0.95_counts_and_l2r.blfix.varfilt.K-30.binarize.FALSE/ldaOut.count_mat.H3K9me3.filt_0.15_0.95_counts_and_l2r.blfix.varfilt.K-30.Robj")
load(inf.lda, v=T)


tm.result <- posterior(out.lda)

dat.umap.k9me3.round2 <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)


# check normal PCA

bad.rname <- "chr8:69800000-69850000"
jmat <- count.mat
jmat <- jmat[rownames(jmat) != bad.rname, ]
jmat.norm <- log2(sweep(jmat, MARGIN = 2, STATS = colSums(jmat), FUN = "/") * 1000 + 1)
jmat.norm2 <- scale(t(jmat.norm), center = TRUE, scale = TRUE)
count.svd <- irlba(A = jmat.norm2, nv = 30, scale = FALSE, center = FALSE)
rownames(count.svd$u) <- rownames(jmat.norm2)
# jmat.norm2.eigengenes <- count.svd$u
dat.svd <- DoUmapAndLouvain(count.svd$u, jsettings = jsettings)
dat.svd.annot <- left_join(dat.svd, subset(dat.metas$H3K9me3, select = c(cell, cluster, clustercol)))

dat.svd.annot2 <- left_join(subset(dat.svd, select = c(cell, umap1, umap2)), subset(dat.metas$H3K9me3, select = c(-umap1, -umap2)))

m.pca <- ggplot(dat.svd.annot, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  ggtitle("PCA") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.pca2 <- ggplot(dat.svd.annot2, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  ggtitle("PCA") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.pca.check <- ggplot(dat.svd.annot2, aes(x = umap1, y = umap2, color = log2(cuts_total))) + 
  geom_point() + 
  scale_color_viridis_c() + 
  ggtitle("PCA") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

dat.umap.k9me3.round2.check <- left_join(dat.umap.k9me3.round2, subset(dat.metas$H3K9me3, select = c(cell, cluster, clustercol, cuts_total)))

m.lda.check <- ggplot(dat.umap.k9me3.round2.check, aes(x = umap1, y = umap2, color = log2(cuts_total))) + 
  geom_point() + 
  scale_color_viridis_c() + 
  ggtitle("LDA") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


lsi.result <- RunLSI(count.mat = as.matrix(jmat), n.components = 30)
dat.lsi.k9me3.round2 <- DoUmapAndLouvain(lsi.result$u, jsettings = jsettings)



outpdf <- file.path(outdir, "lda_vs_pca.H3K9me3_round2only.pdf")
pdf(file = outpdf, useDingbats = FALSE)

m.lda <- ggplot(dat.umap.k9me3.round2 %>% left_join(., dat.metas$H3K9me3 %>% dplyr::select(c(cell, cluster, clustercol))), aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  ggtitle("LDA") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# check raw 



m.lsi <- ggplot(dat.lsi.k9me3.round2 %>% left_join(., dat.metas$H3K9me3 %>% dplyr::select(c(cell, cluster, clustercol))), aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  ggtitle("LSI") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.pca <- ggplot(dat.svd.annot, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  scale_color_identity(guide = "legend", 
                       labels = unique(dat.metas$H3K9me3$cluster)) + 
  ggtitle("PCA") + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

JFuncs::multiplot(m.lda, m.lsi, m.pca, cols = 3)
JFuncs::multiplot(m.lda, m.pca, cols = 2)

lapply(list(m.lda, m.lsi, m.pca), print)

print(m.pca.check)
print(m.lda.check)

JFuncs::multiplot(m.lda.check, m.pca.check, cols = 2)


dev.off()

