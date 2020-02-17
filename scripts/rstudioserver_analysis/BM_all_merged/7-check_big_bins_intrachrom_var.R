# Jake Yeung
# Date of Creation: 2020-01-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/7-check_big_bins_intrachrom_var.R
# Cehck big bins intrachrom var 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)


# Settings ----------------------------------------------------------------

jchromos.auto <- paste("chr", c(seq(19)), sep = "")
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDA  ---------------------------------------------------------------

jmark <- "H3K4me1"
inf <-  paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

tm.result <- posterior(out.lda)
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
dat.var.auto <- CalculateVarAll(dat.impute.log, jchromos.auto)
dat.var <- CalculateVarAll(dat.impute.log, jchromos)


topics.mat <- tm.result$topics
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)


dat.umap.merge <- left_join(dat.umap.long, dat.var.auto)
ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

cells.keep <- dat.umap.merge$cell

# Load raw  ---------------------------------------------------------------


bsize <- "10000000"
if (jmark == "H3K4me1"){
  inf.raw <- paste0("/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_otherWinSize_NoSliding_Blacklist/H3K4me1-BM_SC-merged.tagged.bsize_", bsize, ".step_", bsize, ".countTable.demuxbugfixed.csv")
} else {
  inf.raw <- paste0("/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_otherWinSize_NoSliding_Blacklist/", jmark, "-BM_Linneg_SC-merged.tagged.bsize_", bsize, ".step_", bsize, ".countTable.demuxbugfixed.csv")
}
# inf.raw <- paste0("/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_otherWinSize_NoSliding_Blacklist/", jmark, "-BM_Linneg_SC-merged.tagged.bsize_1000000.step_1000000.countTable.demuxbugfixed.csv")
assertthat::assert_that(file.exists(inf.raw))
mat.raw <- ReadMatSlideWinFormat(inf.raw, as.sparse = TRUE)

mat.raw <- mat.raw[, cells.keep]


# do intrachrom var
mat.norm <- sweep(mat.raw, MARGIN = 2, STATS = colSums(mat.raw), FUN = "/")
# mat.norm <- mat.raw

dat.ncuts <- data.frame(ncuts = apply(mat.raw, 2, sum), cell = colnames(mat.raw), stringsAsFactors = FALSE)

mat.norm.log <- log2(mat.norm * 10^3 + 1)
mat.norm.log.auto <- mat.norm.log[!grepl("^chrX|^chrY", rownames(mat.norm.log)), ]

mat.var <- apply(mat.norm.log.auto, 2, function(jcol) var(jcol, na.rm = TRUE))
dat.var.raw <- data.frame(cell = names(mat.var), var.raw = mat.var, stringsAsFactors = FALSE) %>%
  left_join(dat.ncuts)

# plot(density(mat.var))

dat.merge <- left_join(dat.var.raw, dat.var.auto) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "-"),
         plate = ClipLast(cell, jsep = "_"))

ggplot(dat.merge, aes(y = cell.var.within.sum, x = var.raw, color = experi)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10() + scale_y_log10()  + ggtitle(paste("Bsize:", bsize)) + facet_wrap(~experi)

ggplot(dat.merge, aes(y = var.raw, x = ncuts, color = experi)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10() + scale_y_log10()  + ggtitle(paste("Bsize:", bsize)) + facet_wrap(~experi)

ggplot(dat.merge, aes(y = cell.var.within.sum, x = var.raw)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10() + scale_y_log10() + facet_wrap(~experi) + ggtitle(paste("Bsize:", bsize))


ggplot(dat.merge, aes(y = cell.var.within.sum, x = var.raw, color = experi)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10() + scale_y_log10() + facet_wrap(~experi)

ggplot(dat.merge, aes(y = cell.var.within.sum, x = var.raw, color = plate)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10() + scale_y_log10() + facet_wrap(~plate)

ggplot(dat.merge, aes(x = var.raw, fill = experi, group = experi)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10()
  # scale_x_log10() + scale_y_log10() 
  # facet_wrap(~experi)

ggplot(dat.merge, aes(x = cell.var.within.sum.norm, fill = experi, group = experi)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dat.umap.merge.raw <- left_join(dat.umap.merge, dat.var.raw)

ggplot(dat.umap.merge.raw, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.umap.merge.raw, aes(x = umap1, y = umap2, color = var.raw)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# winsorize?? 

# visualize a high and low var cell?? 


