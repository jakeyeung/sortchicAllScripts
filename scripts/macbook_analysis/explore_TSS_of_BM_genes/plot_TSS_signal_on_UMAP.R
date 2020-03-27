# Jake Yeung
# Date of Creation: 2020-03-17
# File: ~/projects/scchic/scripts/macbook_analysis/explore_TSS_of_BM_genes/plot_TSS_signal_on_UMAP.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


jmark <- "H3K27me3"
jdir <- -1
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load LDA  ---------------------------------------------------------------

inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/LDA_outputs/ldaAnalysisTSS_B6BM_All_allmarks.2020-03-12.var_filt.UnenrichedAndAllMerged.MergedByMarks_final.count_tables_TSS/lda_TSS.", jmark, ".countTableTSS.mapq_40.TSS_50000.blfiltered.K-30.binarize.FALSE/ldaOut.", jmark, ".countTableTSS.mapq_40.TSS_50000.blfiltered.K-30.Robj")
print(inf.lda)
load(inf.lda, v=T)

count.mat.tss <- count.mat

tss.size <- data.frame(cell = colnames(count.mat.tss), tss.size = colSums(count.mat.tss), stringsAsFactors = FALSE)

# Load LDA bins -----------------------------------------------------------

# inf.glm <- "PZ_H3K4me3.AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData"
inf.glm <- paste0("/Users/yeung/data/scchic/from_rstudio/primetime_objs/GLMPCA_outputs.KeepBestPlates2.good_runs/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData")
load(inf.glm, v=T)

bin.size <- data.frame(cell = names(glm.inits$size.factor), bin.size = glm.inits$size.factor, stringsAsFactors = FALSE)

dat.size <- left_join(tss.size, bin.size, by = "cell") %>%
  rowwise() %>%
  mutate(tss.frac = tss.size / bin.size)

# Load annots -------------------------------------------------------------

inf.annot <- paste0("/Users/yeung/data/scchic/from_rstudio/pdfs_all/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf.annot, v=T)

dat.merge <- left_join(dat.umap.glm.fillNAs, dat.size)

pdf(paste0("/Users/yeung/Dropbox/scCHiC_figs/meeting_files/2020-03-17/umap_with_TSS_signal/TSS_signal_on_umap.", jmark, ".pdf"))
  ggplot(dat.merge, aes(x = umap1, y = umap2, color = tss.frac)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = jdir) + ggtitle(jmark, "TSS signal")
  
  ggplot(dat.merge, aes(x = umap1, y = umap2, color = cluster)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark, "TSS signal")
dev.off()
