# Jake Yeung
# Date of Creation: 2020-10-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/4-cleanup_k9me3_mat_rerun_LDA.R
# After figuring out bad clusters, rerun LDA

# rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)


# Load LDA which contains countmat ----------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.allmerged <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld/lda_outputs.count_mat_old_merged_with_new.H3K9me3.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.H3K9me3.K-30.Robj")
load(inf.allmerged, v=T)

# cont.mat.allmerged <- count.mat
count.mat.allmerged <- count.mat

inf.new <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_varfilt.blfix/lda_outputs.count_mat.H3K9me3.filt_0.15_0.95_counts_and_l2r.blfix.varfilt.K-30.binarize.FALSE/ldaOut.count_mat.H3K9me3.filt_0.15_0.95_counts_and_l2r.blfix.varfilt.K-30.Robj"))
load(inf.new, v=T)
count.mat.new <- count.mat


# Load good cells  --------------------------------------------------------

inf.objs <-  file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_scChIX_output/match_UMAP_assign_clusters_objs.2020-10-25.RData")
load(inf.objs, v=T)

ggplot(dat.umap.merged.var.filt.annot %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

good.cells.clean <- subset(dat.umap.merged.var.filt.annot, mark == "H3K9me3" & !is.na(cluster.act))$cell

cells.new.keep <- colnames(count.mat.new) %in% good.cells.clean
print(length(which(cells.new.keep)))
assertthat::assert_that(length(which(cells.new.keep)) > 0)

cells.allmerged.keep <- colnames(count.mat.allmerged) %in% good.cells.clean
print(length(which(cells.allmerged.keep)))
assertthat::assert_that(length(which(cells.allmerged.keep)) > 0)


count.mat.new.filt <- count.mat.new[, cells.new.keep]
count.mat.allmerged.filt <- count.mat.allmerged[, cells.allmerged.keep]

print(dim(count.mat.new.filt))
print(dim(count.mat.allmerged.filt))


outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_H3K9me3_after_filtering")

outf1 <- file.path(outdir, "countmat_H3K9me3_newonly_badclustremoved.rds")
outf2 <- file.path(outdir, "countmat_H3K9me3_allmerged_badclustremoved.rds")

saveRDS(count.mat.new.filt, outf1)
saveRDS(count.mat.allmerged.filt, outf2)


