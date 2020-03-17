# Jake Yeung
# Date of Creation: 2020-02-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_from_hiddendomains/3-MARA_downstream.R
# 
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(hash)
library(igraph)
library(umap)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123



# Load GLMPCA  ------------------------------------------------------------


# load GLMPCA from bins 
jmark <- "H3K27me3"
jexperi <- "AllMerged"
mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1
inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))
load(inf.glm, v=T)

inf.annots <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
assertthat::assert_that(file.exists(inf.annots))
load(inf.annots, v=T)

mdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
assertthat::assert_that(dir.exists(mdir))


outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/mara_outputs"
outname <- paste0("mara_", gsub(".RData", ".pdf", basename(inf.glm)))
outpdf <- file.path(outdir, outname)

pdf(outpdf, useDingbats = FALSE)


dat.glmpca.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)

dat.glmpca.umap <- dat.glmpca.umap %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), "_"))

dat.glmpca.umap$cond <- sapply(dat.glmpca.umap$cell, GetCondFromSamp, mark = jmark)
dat.glmpca.umap$cond <- factor(dat.glmpca.umap$cond, levels = c("Unenriched", "Linneg", "StemCell"))

dat.glmpca.umap.annot <- left_join(dat.glmpca.umap, subset(dat.umap.glm.fillNAs, select = c(cell, cluster, topic.weight)))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m1 <- ggplot(dat.glmpca.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  facet_wrap(~cond)

m2 <- ggplot(dat.glmpca.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  scale_color_manual(values = cbPalette, na.value = "grey85")

print(m1)
print(m2)


# Load MARA output  ------------------------------------------------------

mara.out <- LoadMARA(mdir, make.cnames = FALSE)

# plot top zscores

# integrate motif activity into UMAP and plot 

dat.glmpca.umap.merged <- left_join(dat.glmpca.umap %>% mutate(cell = make.names(cell, unique = FALSE)), mara.out$act.long)

jsub <- subset(mara.out$zscores, zscore > 0.75)

for (i in seq(nrow(jsub))){
  jrow <- jsub[i, ]
  jmotif <- jrow[[1]]
  zscore <- signif(as.numeric(jrow[[2]]), digits = 2)
  print(jmotif)
  m.tmp1 <- ggplot(subset(dat.glmpca.umap.merged, motif == jmotif), aes(x = umap1, y = umap2, color = activity)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    scale_color_viridis_c() + ggtitle(paste(jmotif, "Zscore:", zscore))
  m.tmp2 <- ggplot(subset(dat.glmpca.umap.merged, motif == jmotif), aes(x = umap1, y = umap2, color = activity)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    scale_color_viridis_c() + ggtitle(paste(jmotif, "Zscore:", zscore)) + facet_wrap(~cond)
  print(m.tmp1)
  print(m.tmp2)
}

dev.off()





