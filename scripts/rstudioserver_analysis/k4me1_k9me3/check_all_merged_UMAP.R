# Jake Yeung
# Date of Creation: 2020-10-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/check_all_merged_UMAP.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load  -------------------------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")
# joutsuffix <- "UnionRows_KeepAllCells"; remove.na = FALSE
joutsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters"; remove.na = TRUE
hubprefix <- "/home/jyeung/hub_oudenaarden"

# inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/LDA_outputs/unmix_split_MF_BM_SetupObjs_AllMerged_UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters_H3K4me1-H3K9me3_split_reads/afterUnmixing.MF_SetupObjs_AllMerged_UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters_clstr_by_louvain.H3K4me1xH3K9me3-merged_mat.H3K9me3.K-30.binarize.FALSE/ldaOut.MF_SetupObjs_AllMerged_UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters_clstr_by_louvain.H3K4me1xH3K9me3-merged_mat.H3K9me3.K-30.Robj")
inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/LDA_outputs/unmix_split_MF_BM_SetupObjs_AllMerged_", joutsuffix, "_H3K4me1-H3K9me3_split_reads/afterUnmixing.MF_SetupObjs_AllMerged_", joutsuffix, "_clstr_by_louvain.H3K4me1xH3K9me3-merged_mat.H3K9me3.K-30.binarize.FALSE/ldaOut.MF_SetupObjs_AllMerged_", joutsuffix, "_clstr_by_louvain.H3K4me1xH3K9me3-merged_mat.H3K9me3.K-30.Robj"))

load(inf, v=T)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics

dat.umap <- DoUmapAndLouvain(topics.mat, jsettings)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap$stain <- sapply(dat.umap$cell, function(x) ifelse(grepl("K4me1", x), "dbl", "single"))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = stain)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~stain) + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load fit outputs --------------------------------------------------------


# joutsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters"

# pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_scChIX_output/AllMerged_", joutsuffix, ".pdf")
# pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)

hubprefix <- "/home/jyeung/hub_oudenaarden"

projmain <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/projects_after_unmixing.H3K4me1xH3K9me3/SetupObjs_AllMerged_", joutsuffix))

inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.Robj"
assertthat::assert_that(file.exists(inf.dbl))

# inf.input <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/mouse_spikein_BMround2all.dbl_common_rows_match_dbl/mouse_spikein_BMround2all.dbl_common_rows_match_dbl_clstr_by_louvain_H3K4me1xH3K9me3.removeNA_TRUE.RData"
inf.input <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_", joutsuffix, "/SetupObjs_AllMerged_", joutsuffix, ".clstr_by_louvain_H3K4me1xH3K9me3.removeNA_", remove.na, ".RData")
# inf.output <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/mouse_spikein_BMround2all.dbl_common_rows_match_dbl/unmix_mouse_spikein_BMround2all.dbl_common_rows_match_dbl_clstr_by_louvain_H3K4me1xH3K9me3.removeNA_TRUE.RData"
inf.output <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/SetupObjs_AllMerged_", joutsuffix, "/unmix_SetupObjs_AllMerged_", joutsuffix, ".clstr_by_louvain_H3K4me1xH3K9me3.removeNA_", remove.na, ".RData")


assertthat::assert_that(file.exists(inf.input))
assertthat::assert_that(file.exists(inf.output))

inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_", joutsuffix, "/SetupObjs_AllMerged_", joutsuffix, ".clstr_by_louvain_H3K4me1xH3K9me3.removeNA_", remove.na, ".RData")
# indir.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl"
# assertthat::assert_that(dir.exists(inf.annot))
assertthat::assert_that(file.exists(inf.annot))
# inf.annot.k9me3 <- file.path(indir.annot, "cluster_tables_H3K9me3_BM_all_round2.txt")
# inf.annot.k4me1 <- file.path(indir.annot, "cluster_tables_H3K4me1_BM_all_round2.txt")

load(inf.annot, v=T)

load(inf.input, v=T)
load(inf.output, v=T)




fits.out <- act.repress.coord.lst

w.lst <- sapply(fits.out, function(x) x$w)

# remove 0.01 or 0.99
cells.remove.i <- which(w.lst >= 0.99 | w.lst <= 0.01)
if (length(cells.remove.i) > 0){
  cells.remove <- names(w.lst)[cells.remove.i]
  fits.out[[cells.remove]] <- NULL
}


# if louvains are now from clusters need eto rethink jcoord
cell.vec <- names(fits.out)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- fits.out[[jcell]]
  jweight <- fits.out[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)
  
  # rows are active, columns are repress I THINK?
  # TODO: assumes underscores be careful!
  jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
  jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
  
  if (grepl("_", jlouv.act)){
    jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
  } 
  if (grepl("_", jlouv.repress)){
    jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
  } 
  
  out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()


# coords.dbl.annots <- left_join(coords.dbl, annots.dat)
coords.dbl.annots <- coords.dbl

dat.umap.dbl.merge <- left_join(dat.umap, coords.dbl.annots) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

ggplot(dat.umap.dbl.merge, aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85")

ggplot(dat.umap.dbl.merge, aes(x = umap1, y = umap2, color = exp(lnprob))) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

ggplot(dat.umap.dbl.merge, aes(x = umap1, y = umap2, color = louv.repress)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85")


# Check variance?  --------------------------------------------------------

inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.Robj"
load(inf.dbl, v=T)
tm.result.dbl <- posterior(out.lda)

dat.impute.log <- t(log2(tm.result.dbl$topics %*% tm.result.dbl$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.dbl.merge.var <- left_join(dat.umap.dbl.merge, dat.var)

ggplot(dat.umap.dbl.merge.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

