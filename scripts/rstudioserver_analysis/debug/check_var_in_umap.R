# Jake Yeung
# Date of Creation: 2020-01-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/check_var_in_umap.R
# Try to get out var

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
library(JFuncs)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


#' ## Introduction
#' 
#' The intrachromosomal variance that we are seeing may be coming from differences in MNase efficiency for each well. 
#' If that is the case, we can imagine the UMAP of the data would have different islands representing different celltypes, and then 
#' within each celltype, there would be a gradient of intrachromosomal variance. And there may be many genomic regions 
#' that constitute this additional background signal. But hopefully we can spot a global trend across all the celltypes and then 
#' correct for it using simple linear methods. 




# Load data  --------------------------------------------------------------

# inf <- "/home/jyeung/hpc/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_BM_dblmarks.again2.morememory.fullmerge/lda_outputs.count_mat.H3K4me1-H3K27me3.countcutoff_500.TAcutoff_0.25.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1-H3K27me3.countcutoff_500.TAcutoff_0.25.K-30.Robj"
# inf <- "/home/jyeung/hpc/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_BM_dblmarks.again2.morememory.fullmerge/lda_outputs.count_mat.H3K4me1.countcutoff_500.TAcutoff_0.25.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_500.TAcutoff_0.25.K-30.Robj"
# inf <- "/home/jyeung/hpc/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_BM_dblmarks/lda_outputs.count_mat.BM-B6-H3K4m1-191202.cutoff_500.cutoffTA_0.25.K-30.binarize.FALSE/ldaOut.count_mat.BM-B6-H3K4m1-191202.cutoff_500.cutoffTA_0.25.K-30.Robj"
# inf <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me1.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
# inf <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"

jprefix <- "AllMerged"
jmark <- "H3K4me3"
jdate <- "2019-12-05"

"lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE"

# inf <- "/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_Unenriched_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj"
inf <- paste0("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_", jprefix, "_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.", jdate, ".K-30.binarize.FALSE/ldaOut.B6BM_", jprefix, "_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.", jdate, ".K-30.Robj")

# print(inf)
# print("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.FALSE/ldaOut.B6BM_AllMerged_H3K4me3.TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

inf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
annot.out <- AnnotateBins2(terms.mat = terms.mat, top.thres = 0.995, inf.tss = inf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap.long, dat.var)

#' Here we see that there are 6 islands, each have some gradient of intrachromosomal variance:
ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# plot total counts
dat.counts <- data.frame(cell = colnames(count.mat), cellsize = colSums(count.mat) / 5)
dat.merge <- left_join(dat.merge, dat.counts)

# is column-wise sparsiity of matrix an indicator? 
dat.sparse <- data.frame(cell = colnames(count.mat), frac.zeros = apply(count.mat, 2, function(jcol) nnzero(jcol) / length(jcol)))

dat.merge <- left_join(dat.merge, dat.sparse)

dat.merge <- dat.merge %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         prefix = ClipLast(cell, jsep = "-"))

#' There may be a relationship with the number of cuts, but the relationship isnt as obvious as intrachrom var
ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(cellsize))) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = 1)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = frac.zeros)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = 1)

#' Relationship between cellsize and intrachrom var is weak, but we see Plate 1 is slightly higher than the others. 

ggplot(dat.merge, aes(x = log10(cellsize), y = frac.zeros, color = experi)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge, aes(x = frac.zeros, y = cell.var.within.sum.norm, color = experi)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge, aes(x = log10(cellsize), y = cell.var.within.sum.norm, color = experi)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = frac.zeros, fill = experi)) + geom_density(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~prefix, ncol = 1)

ggplot(dat.merge, aes(x = frac.zeros, fill = experi)) + geom_density(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~prefix, ncol = 1)

ggplot(dat.merge, aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~prefix, ncol = 1)


# plate effect?
m.size <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(cellsize))) + facet_wrap(~experi) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +  
  scale_color_viridis_c(direction = 1)

m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) +  facet_wrap(~experi) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_viridis_c(direction = -1)
  

# ggplot(dat.merge %>% filter(grepl("^PZ", experi)), aes(x = log10(cellsize), y = cell.var.within.sum.norm, color = experi)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#' Lets compare the cellsize and intrachrom var together, separateed by plates. The plate effect seems to be fairly noticeable
print(m.size)
print(m.var)

# multiplot(m.size, m.var, cols = 2)

# can we find PCs that eliminate the cellsize effect?

topics.mat.input <- scale(topics.mat, center = TRUE, scale = TRUE)
scale.fac <- attr(topics.mat.input, "scaled:scale")  # keep my scale value , propagate to terms mat 
const <- attr(topics.mat.input, "scaled:center")
const.terms <- sweep(terms.mat, MARGIN = 1, STATS = const, FUN = "*")
terms.mat.input <- sweep(terms.mat, MARGIN = 1, STATS = scale.fac, FUN = "*")
impute.input <- t(topics.mat.input %*% terms.mat.input)  # currently I am off by a constant...
# back out the constant
# impute.input <- impute.input + const.terms

dat.comp <- prcomp(topics.mat.input, center = FALSE, scale. = FALSE)

dat.pca <- data.frame(cell = rownames(dat.comp$x), dat.comp$x, stringsAsFactors = FALSE)

dat.merge.pca <- left_join(dat.merge, dat.pca)

m.louv.pca <- ggplot(dat.merge.pca, aes(x = PC1, y = PC2, color = louvain)) + geom_point() + scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# m.louv.pca <- ggplot(dat.merge.pca, aes(x = PC2, y = PC3, color = louvain)) + geom_point() + scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Usually we look at the cell-to-topics matrix by visualizing on UMAP. let's look at the cell-to-topics matrix in PCA space instead...

print(m.louv)
print(m.louv.pca)

#' first PC maybe shows the intrachromosomal variability more cleanly:

ggplot(dat.merge.pca, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge.pca, aes(x = PC1, y = PC2, color = log10(cellsize))) + geom_point() + scale_color_viridis_c(direction = 1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#' Higher order PCs probably related to celltypes without the intrachromosomal variability...
ggplot(dat.merge.pca, aes(x = PC2, y = PC3, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge.pca, aes(x = PC2, y = PC3, color = louvain)) + geom_point() + scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.merge.pca, aes(x = PC3, y = PC4, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge.pca, aes(x = PC3, y = PC4, color = louvain)) + geom_point() + scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.merge.pca, aes(x = PC4, y = PC5, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge.pca, aes(x = PC4, y = PC5, color = louvain)) + geom_point() + scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#' First PC seems to be related to intrachromosomal variability... pretty strong effect and it's all on PC1, this is promising
ggplot(dat.merge.pca, aes(x = PC1, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#' higher order PCs not strongly related to intrachrom var, that's good
ggplot(dat.merge.pca, aes(x = PC2, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge.pca, aes(x = PC3, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge.pca, aes(x = PC4, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#' Let's remove PC1 and keep components 2 to 15...
plot(dat.comp$sdev ^ 2 / sum(dat.comp$sdev), xlab = "PC", ylab = "Frac of Var Explained")  # 30 percent of variance coming from PC1!!

rot.mat <- dat.comp$rotation[, 2:15]  # use this for genes later
# rot.mat <- dat.comp$rotation  # use this for genes later

dat.proj <- data.frame(topics.mat.input %*% rot.mat, stringsAsFactors = FALSE)

#' We now have a new gene to topics matrix, where we removed a principal component that contained a strong effect we think may be an artifact:
#' let's put this new gene-to-topics matrix into the UMAP and see what it looks like now

umap.out <- umap(dat.proj, config = jsettings)
dat.umap.proj <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)

# add meta dat
dat.umap.proj <- left_join(dat.umap.proj, subset(dat.merge, select = c(-umap1, -umap2)))

#' We see now intrachrom var is not a strong effect in the UMAP
ggplot(dat.umap.proj, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle("UMAP after removing PC1")

m.umap.corrected <- ggplot(dat.umap.proj, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

#' Lets compare the same louvain clustering before and after

print(m.louv + ggtitle("UMAP before correction"))
print(m.umap.corrected + ggtitle("UMAP after correction. Louvain matches original UMAP"))

# multiplot(m.louv, m.umap.corrected)

# check plate effect
m.umap.corrected.plate <- ggplot(dat.umap.proj, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)

#' The plates differences may have improved, but not completely solved...

print(m.var + facet_wrap(~experi) + ggtitle("Plate effect before correction"))
print(m.umap.corrected.plate + ggtitle("Plate effect after correction"))

dat.sum <- dat.merge %>%
  group_by(louvain, experi) %>%
  summarise(ncell = length(cell))

m.counts.louv.old <- ggplot(dat.sum, aes(x = louvain, group = experi, fill = experi, y = ncell)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 

# redo louvain in new space

dat.proj.louvain <- DoLouvain(dat.proj, jsettings, subset(dat.umap.proj, select = -louvain))

m.louv.new <- ggplot(dat.proj.louvain, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

# multiplot(m.louv, m.umap.corrected)

print(m.louv + ggtitle("UMAP before correction"))
print(m.louv.new + ggtitle("UMAP after correction. Louvain was recalculated on corrected matrix"))


# check counts in new space
dat.sum.new <- dat.proj.louvain %>%
  group_by(louvain, experi) %>%
  summarise(ncell = length(cell))

m.counts.louv.new <- ggplot(dat.sum.new, aes(x = louvain, group = experi, fill = experi, y = ncell)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 

multiplot(m.counts.louv.old, m.counts.louv.new)

# ggplot(subset(dat.merge, louvain %in% c(4, 6)), aes(x = umap1, y = umap2, color = louvain)) + geom_point()   + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_grid(louvain ~ experi)
# 
# ggplot(subset(dat.proj.louvain, louvain %in% c(6, 2)), aes(x = umap1, y = umap2, color = louvain)) + geom_point()   + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_grid(louvain ~ experi)

# Plot expression of a gene onto this umap by propagating the rotation matrix 
dim(tm.result$terms)
# dim(tm.result$topics)

genes.proj <- t(terms.mat.input) %*% rot.mat
# check sums after a truncated projection
colSums(genes.proj)  # sums to 0

# check loadings that we removed

loadings.var <- data.frame(term = colnames(terms.mat.input), loading = t(terms.mat.input) %*% dat.comp$rotation[, 1], stringsAsFactors = FALSE) %>%
  ungroup() %>%
  mutate(rnk = rank(-loading)) %>%
  arrange(rnk)

#' Explore the suspicious peaks... 
ggplot(loadings.var %>% filter(rnk < 50), aes(x = rnk, y = log10(loading), label = term)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#' As expected, suspicious peaks have large loadings for topics that are also in PC1
sort(log10(terms.mat.input[, "chr17:35160000-35260000"]), decreasing = TRUE)
sort(dat.comp$rotation[, 1], decreasing = TRUE)  

# plot a gene onto the umap
dat.impute <- genes.proj %*% t(dat.proj)

# plot S100a8, granulocytes
jgenes <- c("Ccl5", "Sox6", "Bach2", "Ccl2", "S100a8")

jterms.dat <- subset(annot.out$out2.df.closest, gene %in% jgenes) %>%
  group_by(gene) %>%
  filter(abs(dist.to.tss) == min(abs(dist.to.tss)))
jterms <- jterms.dat$region_coord
jgenes.ordered <- jterms.dat$gene
jgenes.hash <- hash(make.names(jterms), jgenes.ordered)

exprs <- data.frame(cell = colnames(dat.impute), t(dat.impute[jterms, ]), stringsAsFactors = FALSE)
# rename columns?
colnames(exprs) <- sapply(colnames(exprs), function(x) ifelse(!is.null(jgenes.hash[[x]]), jgenes.hash[[x]], x))

dat.withgene <- left_join(dat.umap.proj, exprs)

#' Lets project the genes using the same rotation matrix and then create a corrected expression matrix
for (jgene in jgenes){
  print(jgene)
  m <- PlotXYWithColor(dat.withgene, xvar = "umap1", yvar = "umap2", cname = jgene, jtitle = paste("Imputed exprs of", jgene))
  print(m)
}




# # recalculate intra chrom var?
# dat.var.proj <- CalculateVarAll(log2(dat.impute), jchromos)
# dat.merge.proj <- left_join(dat.proj, dat.var.proj)

# plot onto new umap



# # are these from toipc_1? Not enough, the granu trajectory has a large var-effect probably
# umap.out.notopic1 <- umap(topics.mat[, -1], config = jsettings)
# dat.umap.long.notopic1 <- data.frame(cell = rownames(umap.out.notopic1[["layout"]]), 
#                                      umap1 = umap.out.notopic1[["layout"]][, 1], 
#                                      umap2 = umap.out.notopic1[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long.notopic1 <- left_join(dat.umap.long.notopic1, dat.merge %>% select(., c(-umap1, -umap2)))
# ggplot(dat.umap.long.notopic1, aes(x = umap1, y = umap2, color = louvain)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette)
# 



