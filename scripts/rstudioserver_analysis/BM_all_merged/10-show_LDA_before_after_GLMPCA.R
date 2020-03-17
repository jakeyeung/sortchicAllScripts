# Jake Yeung
# Date of Creation: 2020-02-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/10-show_LDA_before_after_GLMPCA.R
# Compare GLMPCA before and after

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)
library(topicmodels)

library(scchicFuncs)


# Load data ---------------------------------------------------------------

jmark <- "H3K9me3"
jsuffix <- "Unenriched"
# inf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs/PZ_H3K4me1.KeepMorePlates.GLMPCA_var_correction.150.2020-02-02.binskeep_1000.devmode.RData"
# inf <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs/PZ_", jmark, ".GLMPCA_var_correction.150.2020-02-01.binskeep_500.devmode.RData")
# inf <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs/PZ_", jmark, ".GLMPCA_var_correction.150.2020-02-01.binskeep_500.devmode.RData")
inf <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs/PZ_", jmark, ".Unenriched.GLMPCA_var_correction.150.2020-02-04.binskeep_1000.devmode..RData")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

unique(sapply(colnames(glm.inits$Y.filt), function(x) ClipLast(x, jsep = "_")))

# inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt/lda_outputs.BM_", jmark, ".varcutoff_0.3.platesRemoved.K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.platesRemoved.K-30.Robj")
# inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt/lda_outputs.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.K-30.Robj")
inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt/lda_outputs.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.", jsuffix, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.", jsuffix, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.lda))

load(inf.lda, v=T)
 

# Show UMAP before correction ---------------------------------------------

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long) %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), jsep = "_"))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = plate)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = plate)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~plate)


dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~plate)

# GLM after
topics.mat.glm <- glm.out$factors

umap.out <- umap(topics.mat.glm, config = jsettings)
dat.umap.long.glm <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE) %>% 
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), jsep = "_")) %>%
  left_join(., dat.var) %>%
  left_join(., dat.umap.long %>% dplyr::select(cell, louvain))

ggplot(dat.umap.long.glm, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~plate)

m.before <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + ggtitle(jmark)
m.after <- ggplot(dat.umap.long.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) 

multiplot(m.before, m.after)



# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)










