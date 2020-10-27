# Jake Yeung
# Date of Creation: 2020-10-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/5-downstream_connect_UMAPs_together.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#851663", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load data ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inobjs <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_scChIX_output/match_UMAP_assign_clusters_objs.2020-10-25.RData")

load(inobjs, v=T)


# Plot nice umaps ---------------------------------------------------------

dat.umap.merged.var.filt.annot.grey <- dat.umap.merged.var.filt.annot %>%
  rowwise() %>%
  mutate(cluster.act = ifelse(stain == "single", NA, cluster.act)) %>%
  group_by(mark) %>%
  mutate(umap2.scale = scale(umap2, center = TRUE, scale = TRUE))

ggplot(dat.umap.merged.var.filt.annot.grey %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point()  + 
  facet_wrap(~stain) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged.var.filt.annot.grey, aes(x = umap1.shift, y = umap2.scale, color = cluster.act, group = cell)) + 
  geom_point()  + 
  geom_line(alpha = 0.05) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.umap.merged.var.filt.annot %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster.act, group = cell)) + 
  geom_point()  + 
  geom_line(alpha = 0.05) + 
  facet_wrap(~stain) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# how many double cells
length(unique(subset(dat.umap.merged.var.filt.annot.grey, stain == "dbl")$cell))

# Check after removing ----------------------------------------------------

inf.removed <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_badclustremoved/lda_outputs.countmat_H3K9me3_newonly_badclustremoved.K-30.binarize.FALSE/ldaOut.countmat_H3K9me3_newonly_badclustremoved.K-30.Robj")
load(inf.removed, v=T)

tm.result.rm <- posterior(out.lda)

dat.umap.rm <- DoUmapAndLouvain(tm.result.rm$topics, jsettings)

louv2clst <- list("1" = "HSPCs", 
                  "2" = "Myeloid", 
                  "3" = "Myeloid",
                  "4" = "Eryth",
                  "5" = "Lymphoid")

dat.umap.rm$cluster <- sapply(dat.umap.rm$louvain, function(x) AssignHash(x, louv2clst, null.fill = NA))

dat.umap.rm$plate <- sapply(dat.umap.rm$cell, function(x) ClipLast(x, jsep = "_"))
dat.umap.rm$experi <- sapply(dat.umap.rm$plate, function(x) ClipLast(x, jsep = "-"))

ggplot(dat.umap.rm, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  ggtitle("H3K9me3 cleaned") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.rm, aes(x = umap1, y = umap2, color = plate)) + 
  geom_point(alpha = 0.25) + 
  ggtitle("H3K9me3 cleaned") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.raw.allmerged <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_H3K9me3_after_filtering/countmat_H3K9me3_allmerged_badclustremoved.rds")

mat.allmerged <- readRDS(inf.raw.allmerged)

cuts.dat <- data.frame(cell = colnames(mat.allmerged), ncuts = colSums(mat.allmerged), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(plate = ifelse(grepl("rep", cell), "New", "Old"))


ggplot(cuts.dat, aes(x = log10(ncuts), fill = plate)) + 
         geom_density(alpha = 0.25) + 
  theme_bw() + 
  ggtitle("Cuts between two plates") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check after removing AllMerged ------------------------------------------



inf.removed2 <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_badclustremoved/lda_outputs.countmat_H3K9me3_allmerged_badclustremoved.K-30.binarize.FALSE/ldaOut.countmat_H3K9me3_allmerged_badclustremoved.K-30.Robj")
load(inf.removed2, v=T)

tm.result.rm.allmerged <- posterior(out.lda)

dat.umap.rm.allmerged <- DoUmapAndLouvain(tm.result.rm.allmerged$topics, jsettings)

dat.umap.rm.allmerged$plate <- sapply(dat.umap.rm.allmerged$cell, function(x) ClipLast(x, jsep = "_"))
dat.umap.rm.allmerged$experi <- sapply(dat.umap.rm.allmerged$plate, function(x) ClipLast(x, jsep = "-"))

ggplot(dat.umap.rm.allmerged, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  ggtitle("H3K9me3 cleaned AllMerged") + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.rm.allmerged, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  facet_wrap(~experi) + 
  ggtitle("H3K9me3 cleaned AllMerged") + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load annots -------------------------------------------------------------

dat.annot.k9me3 <- subset(dat.umap.merged.var.filt.annot, mark == "H3K9me3" & stain == "single" & !is.na(cluster.act), select = c(cell, cluster.act))

dat.umap.rm.allmerged.annot <- left_join(dat.umap.rm.allmerged, dat.annot.k9me3)

louv2clst.allmerged <- list("1" = "HSPCs", 
                            "2" = "Myeloid", 
                            "3" = "Eryth",
                            "4" = "Lymphoid",
                            "5" = "HSPCs",
                            "6" = "Myeloid")

dat.umap.rm.allmerged.annot$cluster <- sapply(dat.umap.rm.allmerged.annot$louvain, function(x) AssignHash(x, louv2clst.allmerged, null.fill = NA))


ggplot(dat.umap.rm.allmerged.annot, aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  ggtitle("H3K9me3 cleaned AllMerged") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.rm.allmerged.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(alpha = 0.8) + 
  ggtitle("H3K9me3 cleaned All Merged") + 
  facet_wrap(~experi) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.rm.allmerged.annot, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point(alpha = 0.8) + 
  ggtitle("H3K9me3 cleaned All Merged") + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

