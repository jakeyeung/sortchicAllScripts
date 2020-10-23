# Jake Yeung
# Date of Creation: 2020-08-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse_LDA_downstream.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


hubprefix <- "/home/jyeung/hub_oudenaarden"


inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein/lda_outputs.H3K4me3_merged_with_old.K-30.binarize.FALSE/ldaOut.H3K4me3_merged_with_old.K-30.Robj")

load(inf, v=T)
count.mat.new <- count.mat; rm(count.mat)
out.lda.new <- out.lda

tm.result <- posterior(out.lda.new)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(posterior(out.lda.new)$topics, jsettings = jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))

jchromos <- paste("chr", seq(19), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.merge <- left_join(dat.umap, dat.var)

ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)


inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.merged_old_and_new/cell_cluster_merged_with_spikein_info.txt")

dat.annot <- fread(inf.annot)

dat.umap.merge2 <- left_join(subset(dat.umap.merge, select = -c(louvain)), subset(dat.annot, select = -c(umap1, umap2)))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.merge2, aes(x = umap1, y = umap2, color = cluster.infer)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cond)

ggplot(dat.umap.merge2, aes(x = umap1, y = umap2, color = )) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cond)

ReadMatTSSFormat(inf = , as.sparse = , add.coord = , sort.rnames = )

# 
# 
# # Load existing manifold --------------------------------------------------
# 
# 
# inf.old <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
# 
# load(inf.old, v=T)
# 
# rownames.keep <- rownames(count.mat)
# rownames.keep.i <- which(rownames(count.mat.new) %in% rownames.keep)
# rownames.keep.innew <- rownames(count.mat.new)[rownames.keep.i]
# rownames.keep.toadd <- rownames(count.mat)[which(!rownames(count.mat) %in% rownames.keep.innew)]
# 
# 
# mat.filt <- count.mat.new[rownames.keep.i, ]
# # add zeros to nonexisting mats
# mat.zeros <- matrix(data = 0, ncol = ncol(mat.filt), nrow = length(rownames.keep.toadd), dimnames = list(rownames.keep.toadd, colnames(mat.filt)))
# 
# mat.filt2 <- Matrix(rbind(mat.filt, mat.zeros)[rownames.keep, ], sparse = TRUE)
# 
# assertthat::assert_that(identical(rownames(count.mat), rownames(mat.filt2)))
# 
# # write eto output
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.for_projections"
# saveRDS(mat.filt2, file = file.path(outdir, paste0("H3K4me3_padded_zeros_for_projections.rds")))
# 
# # merge with previous mat for LDA 
# 
# mat.merged <- cbind(count.mat, mat.filt2)
# saveRDS(mat.merged, file = file.path(outdir, paste0("H3K4me3_merged_with_old.rds")))
# 
# 
# 





