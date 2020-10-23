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


inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein/lda_outputs.H3K4me3_BM.K-30.binarize.FALSE/ldaOut.H3K4me3_BM.K-30.Robj")
# inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_ze")
# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_zebrafish_spikein/lda_outputs.H3K4me3_WKM.K-30.binarize.FALSE/ldaOut.H3K4me3_WKM.K-30.Robj"

load(inf, v=T)
count.mat.new <- count.mat; rm(count.mat)
out.lda.new <- out.lda


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(posterior(out.lda.new)$topics, jsettings = jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load existing manifold --------------------------------------------------


inf.old <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")

load(inf.old, v=T)

rownames.keep <- rownames(count.mat)
rownames.keep.i <- which(rownames(count.mat.new) %in% rownames.keep)
rownames.keep.innew <- rownames(count.mat.new)[rownames.keep.i]
rownames.keep.toadd <- rownames(count.mat)[which(!rownames(count.mat) %in% rownames.keep.innew)]


mat.filt <- count.mat.new[rownames.keep.i, ]
# add zeros to nonexisting mats
mat.zeros <- matrix(data = 0, ncol = ncol(mat.filt), nrow = length(rownames.keep.toadd), dimnames = list(rownames.keep.toadd, colnames(mat.filt)))

mat.filt2 <- Matrix(rbind(mat.filt, mat.zeros)[rownames.keep, ], sparse = TRUE)

assertthat::assert_that(identical(rownames(count.mat), rownames(mat.filt2)))

# write eto output
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.for_projections"
saveRDS(mat.filt2, file = file.path(outdir, paste0("H3K4me3_padded_zeros_for_projections.rds")))

# merge with previous mat for LDA 

mat.merged <- cbind(count.mat, mat.filt2)
saveRDS(mat.merged, file = file.path(outdir, paste0("H3K4me3_merged_with_old.rds")))








