# Jake Yeung
# Date of Creation: 2020-01-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/4-explore_bone_marrow_progenitors.R
# Find these darn granulocyte progenitors, maybe integrate with Giladi et al. 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

# Load data  --------------------------------------------------------------

# inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"
# inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"
inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16/lda_outputs.B6BM_Unenriched_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.K-50.Robj"
assertthat::assert_that(file.exists(inf))
load(inf, v=T)


inf.raw <- "/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_TSS/H3K4me3-BM_Linneg_SC-merged.tagged.countTable.TSS_50000.csv"

count.tss <- ReadMatTSSFormat(inf.raw, as.sparse = TRUE)
count.tss.long <- CollapseRowsByGene(count.tss, as.long = TRUE)
colnames(count.tss.long) <- c("gene", "cell", "count")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

topics.mat <- posterior(out.lda)$topics

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# add plate and experi
dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "-"),
         plate = ClipLast(cell, jsep = "_"))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette)

print(table(dat.umap.long$experi))

# ggplot(dat.umap.long %>% filter(experi == "B6-13W1-BM-H3K4me3"), aes(x = umap1, y = umap2, color = plate)) + 
#   geom_point() + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
#   scale_color_manual(values = cbPalette)


# ggplot(dat.umap.long %>% filter(experi == "PZ-Bl6-BM-Linneg-H3K4me3"), aes(x = umap1, y = umap2, color = plate)) + 
#   geom_point() + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
#   scale_color_manual(values = cbPalette)

# find Elane??

jgenes <- c("Elane", "F13a1", "Irf8", "Cebpe", "S100a9", "S100a8", "Hlf")

# expression of genes
count.tss

jsub <- data.table::dcast(data = count.tss.long %>% group_by(cell) %>% mutate(count = count / sum(count)), 
                          formula = cell ~ gene, value.var = "count")
# as.data.frame(jsub) <- jsub
# rownames(jsub) <- jsub$cell
# jsub$cell <- NULL

dat.umap.long.genes <- left_join(dat.umap.long, jsub)

PlotXYWithColor(dat.umap.long.genes %>% filter(umap2 > 5), xvar = "umap1", yvar = "umap2", cname = "Elane", remove.axis.info = FALSE)

PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "S100a8")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Hlf")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Irf8")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "S100a9")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Ccl5")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Sox6")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Hoxd4")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Cebpb")
PlotXYWithColor(dat.umap.long.genes, xvar = "umap1", yvar = "umap2", cname = "Cebpe")


# signal noise ratio for plates -------------------------------------------

count.sum <- data.frame(cell = colnames(count.mat), counts = colSums(count.mat), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(cell, jsep = "-"))

print(table(count.sum$experi))

ggplot(count.sum %>% filter(experi == "B6-13W1-BM-H3K4me3"), aes(x = log10(counts), fill = plate)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(count.sum %>% filter(experi == "B6-13W1-BM-H3K4me3"), aes(x = log10(counts), fill = plate)) + geom_histogram() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~plate, ncol = 1)

# subset granulocytes
cells.sub <- unique(subset(dat.umap.long.genes, umap2 > 5)$cell)

ggplot(count.sum %>% filter(experi == "B6-13W1-BM-H3K4me3" & cell %in% cells.sub), aes(x = log10(counts), fill = plate)) + geom_histogram() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~plate, ncol = 1)



dat.merge <- left_join(subset(dat.umap.long.genes, select = c(cell, umap1, umap2, louvain, experi, plate, Elane)), subset(count.sum, select = c(cell, counts)))

PlotXYWithColor(subset(dat.merge, experi == "B6-13W1-BM-H3K4me3" & cell %in% cells.sub), xvar = "umap1", yvar = "umap2", cname = "counts", remove.axis.info = FALSE)

ggplot(dat.merge, aes(x = counts, y = Elane)) + geom_point()

count.mat.sub <- count.mat[, cells.sub]
bin.sum <- data.frame(cell = rownames(count.mat.sub), counts = rowSums(count.mat.sub), stringsAsFactors = FALSE) %>%
  rowwise()



