# Jake Yeung
# Date of Creation: 2020-01-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-analyze_LDA_downstream_with_celltypes_gibbs.correct_background.R
# Correct background 

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

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(JFuncs)
library(scchicFuncs)


# Load data ---------------------------------------------------------------


# inf <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.NoScraping.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.binarize.FALSE/ldaOut.mat.NoScraping.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.Robj"
inf <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.Robj"
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

topics.mat <- posterior(out.lda)$topics
terms.mat <- posterior(out.lda)$terms

annots.out <- AnnotateBins(terms.mat, top.thres = 0.995, inf.tss = "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

dat.impute.log <- log2(t(topics.mat %*% terms.mat))
dat.impute.lin <- t(topics.mat %*% terms.mat)
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

# Regress out  ------------------------------------------------------------

topics.mat.scaled <- scale(topics.mat, center = TRUE, scale = TRUE)
terms.mat.scaled <- scale(terms.mat, center = TRUE, scale = TRUE)

topics.pca <- prcomp(topics.mat.scaled, center = FALSE, scale. = FALSE)

plot(topics.pca$sdev ^ 2 / sum(topics.pca$sdev ^ 2))

dat.pca <- data.frame(cell = rownames(topics.mat), data.frame(topics.pca$x), stringsAsFactors = FALSE)

dat.merge <- left_join(dat.umap.long, dat.var)
dat.merge <- left_join(dat.merge, dat.pca)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.merge, aes(x = PC1, y = PC2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.merge, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)
 
ggplot(dat.merge, aes(x = -PC1, y = cell.var.within.sum.norm, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.merge, aes(x = -PC2, y = cell.var.within.sum.norm, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)



# From PC1 ----------------------------------------------------------------

screeplot(topics.pca)
plot(topics.pca$sdev ^ 2 / sum(topics.pca$sdev ^ 2))
rot.mat <- topics.pca$rotation[, 2:15]

topics.mat.proj <- topics.mat.scaled %*% rot.mat
terms.mat.proj <- t(terms.mat.scaled) %*% rot.mat

# show corrected umap

umap.out.proj <- umap(topics.mat.proj, config = jsettings)
dat.umap.long.proj <- data.frame(cell = rownames(umap.out.proj[["layout"]]), 
                                 umap1 = umap.out.proj[["layout"]][, 1], 
                                 umap2 = umap.out.proj[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long.proj <- DoLouvain(topics.mat.proj, jsettings, dat.umap.long.proj, clstr.cname = "louvain_AfterCorrection")

dat.merge.proj <- left_join(dat.umap.long.proj, dat.var)

# add old louvain
dat.merge.proj <- left_join(dat.merge.proj, subset(dat.umap.long, select = c(cell, louvain)))


m.louv.before <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
m.louv.after <- ggplot(dat.merge.proj, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)


JFuncs::multiplot(m.louv.before, m.louv.after)


# redo louvain
m.louv.after.redo <- ggplot(dat.merge.proj, aes(x = umap1, y = umap2, color = louvain_AfterCorrection)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
print(m.louv.after)


# Plot genes ?  -----------------------------------------------------------

jgenes <- c("Lgr5", "Alpi", "Lrig1", "Dclk1", "Rfx6", "Sox6", "Ikzf1", "Foxo1", "Gimap6", "Bcl11a", "S100a8", "Ebf1", "Elane")

imput.mat <- topics.mat.proj %*% t(terms.mat.proj)


jterms <- subset(annots.out$terms.filt, gene %in% jgenes) %>%
  group_by(gene) %>%
  filter(rnk == min(rnk))
jgenes.terms <- make.names(jterms$termgene)

jsub <- as.data.frame(imput.mat[, jterms$term])
colnames(jsub) <- jgenes.terms

dat.exprs.proj <- data.frame(cell = rownames(imput.mat), jsub, stringsAsFactors = FALSE)
dat.merge.proj.exprs <- left_join(dat.merge.proj, dat.exprs.proj)


pdf("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs/correct_var_LDA/k4me1_LDA_downstream_var_corrected.pdf")

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.merge, aes(x = PC1, y = PC2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.merge, aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.merge, aes(x = -PC1, y = cell.var.within.sum.norm, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.merge, aes(x = -PC2, y = cell.var.within.sum.norm, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

print(m.louv.before)
print(m.louv.after)
print(m.louv.after.redo)

for (jgene in jgenes.terms){
  print(jgene)
  m <- PlotXYWithColor(dat.merge.proj.exprs, xvar = "umap1", yvar = "umap2", cname = jgene)
  print(m)
}

dev.off()

# Check Alpi on original UMAP? 

# imput.mat.orig <- t(dat.impute.log)
imput.mat.orig <- t(dat.impute.lin)
jsub.orig <- as.data.frame(imput.mat.orig[, jterms$term])
colnames(jsub.orig) <- jgenes.terms

dat.exprs.orig <- data.frame(cell = rownames(imput.mat.orig), jsub.orig, stringsAsFactors = FALSE)
dat.merge.orig.exprs <- left_join(dat.merge, dat.exprs.orig)

for (jgene in jgenes.terms){
  print(jgene)
  m <- PlotXYWithColor(dat.merge.orig.exprs, xvar = "umap1", yvar = "umap2", cname = jgene)
  print(m)
}

# plot a topic 
# dat.merge.topics <- left_join(data.frame(cell = rownames(topics.mat), log(topics.mat), stringsAsFactors = FALSE), dat.umap.long)
dat.merge.topics <- left_join(data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE), dat.umap.long)

PlotXYWithColor(dat.merge.topics, xvar = "umap1", yvar = "umap2", cname = "X5")

jterms.dat <- subset(annots.out$terms.filt, topic == 5 & gene == "Dclk1") 
jterm <- jterms.dat$term[[1]]

# exprs.dat <- data.frame(cell = colnames(dat.impute.log), exprs = dat.impute.log[jterm, ], stringsAsFactors = FALSE)
exprs.dat <- data.frame(cell = colnames(dat.impute.lin), exprs = dat.impute.lin[jterm, ], stringsAsFactors = FALSE)

exprs.merge <- left_join(dat.merge.topics, exprs.dat)

PlotXYWithColor(exprs.merge, xvar = "umap1", yvar = "umap2", cname = "exprs")


# rownamesd(dat.blood.long) <- dat.blood$`Gene Name`
# dat.blood.long <- reshape2::melt(dat.blood.long, variable.name = "celltype", value.name = "logexprs")



# Save output -------------------------------------------------------------

outf <- file.path("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/rdata", "k4me1_LDA_downstream_var_correted.RData")
save(outf, topics.mat, terms.mat, rot.mat, topics.mat.proj, terms.mat.proj, annots.out, dat.merge, dat.merge.proj, file = outf)


# write bad cells to output
cells.bad.outf <- file.path("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/cells_list", "k4me1_bad_cells_maybe_immune.txt")
cells.good.outf <- file.path("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/cells_list", "k4me1_good_cells_removed_immune.txt")

ggplot(dat.merge.proj, aes(x = umap1, y = umap2, color = louvain_AfterCorrection)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

cells.remove <- subset(dat.merge.proj, umap2 > 5)$cell
cells.cellskeep <- subset(dat.merge.proj, umap2 < 5)$cell

fwrite(as.data.frame(cells.remove), file = cells.bad.outf, col.names = FALSE)
fwrite(as.data.frame(cells.cellskeep), file = cells.good.outf, col.names = FALSE)




