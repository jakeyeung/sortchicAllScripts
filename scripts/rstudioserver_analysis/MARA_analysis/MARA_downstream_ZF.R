# Jake Yeung
# Date of Creation: 2020-08-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/MARA_downstream_ZF.R
# Analyze MARA downstream for ZF 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load mara outputs -------------------------------------------------------

jmark <- "H3K4me1"
# mdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/countmat_PZ_fromHiddenDomains_H3K4me1.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
mdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/", jmark, "/mara_output/ldaOut.H3K4me1.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.keepNbins_250-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/ldaOut.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.keepNbins_250")
assertthat::assert_that(dir.exists(mdir))

mara.out <- LoadMARA(mdir = mdir, make.cnames = FALSE)

zscores.sub <- subset(mara.out$zscores, zscore > 0.7)
motifs.keep <- zscores.sub$motif


# Load ZF UMAP  -----------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap.lst <- lapply(jmarks, function(jmark){
  inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks/lda_outputs.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.binarize.FALSE/ldaOut.", jmark, ".imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.Robj")
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  
  m1 <- ggplot(dat.umap, aes(x = umap1, y = umap2)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
  print(m1)
  dat.umap$mark <- jmark
  return(dat.umap)
})


# Load annots -------------------------------------------------------------

dat.annot.lst <- lapply(jmarks, function(jmark){
  inf.annot <- file.path(hubprefix, paste0("jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt"))
  assertthat::assert_that(file.exists(inf.annot))
  dat.annot <- fread(inf.annot)
  dat.annot$mark <- jmark
  return(dat.annot)
})

dat.merge.lst <- lapply(jmarks, function(jmark){
  dat.merge <- left_join(dat.umap.lst[[jmark]], subset(dat.annot.lst[[jmark]], select = c(cell, cluster)))
  return(dat.merge)
})

m.lst <- lapply(jmarks, function(jmark){
  dat.merge <- dat.merge.lst[[jmark]]
  m1 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~mark)
  return(m1)
})


# Overlay activities onto UMAP  -------------------------------------------

act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
colnames(act.mat.clean) <- mara.out$act.mat$motif
act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
  ungroup() %>%
  mutate(cell = gsub("\\.", "-", cell))


jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
rownames(jmat) <- act.mat.clean.dat$cell

jmat.sub <- jmat[, motifs.keep]

# overlay 

print(motifs.keep)
jmotif <- "Ikzf1"
jmotif <- "Erg"
jmotif <- "Spic"
jmotif <- "Spib"
jmotif <- "Nkx3.2"
jmotif <- "Hlf"
jmotif <- "Cebpb"
jmotif <- "Sox9"
jmotif <- "Irf4"
jmotif <- "Mef2c"


jmotif <- "Etv4"

jmotif <- "Runx1"
jmotif <- "Gata1"
jmotif <- "Tal1"
jmotif <- "Hoxa1"

jsub.dat <- data.frame(activity = jmat.sub[, jmotif], cell = rownames(jmat.sub), stringsAsFactors = FALSE) %>%
  left_join(., dat.merge.lst[[jmark]])

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(jsub.dat, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(jsub.dat, aes(x = umap1, y = umap2, color = activity)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()




# Check bone marrow  ------------------------------------------------------


mdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/countmat_PZ_fromHiddenDomains_H3K4me1.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
assertthat::assert_that(dir.exists(mdir))

mara.out <- LoadMARA(mdir = mdir, make.cnames = FALSE)

act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
colnames(act.mat.clean) <- mara.out$act.mat$motif
act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
  ungroup() %>%
  mutate(cell = gsub("\\.", "-", cell))

dat.merge <- left_join(dat.glmpca.umap.annot, act.mat.clean.dat)


