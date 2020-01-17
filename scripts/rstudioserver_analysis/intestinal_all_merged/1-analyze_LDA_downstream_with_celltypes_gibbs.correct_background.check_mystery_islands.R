# Jake Yeung
# Date of Creation: 2020-01-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-analyze_LDA_downstream_with_celltypes_gibbs.correct_background.check_mystery_islands.R
# 

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


# Plot genes ?  -----------------------------------------------------------

jgenes <- c("Lgr5", "Alpi", "Lrig1", "Dclk1", "Rfx6", "Sox6", "Ikzf1", "Foxo1", "Gimap6", "Bcl11a")

imput.mat <- topics.mat.proj %*% t(terms.mat.proj)


jterms <- subset(annots.out$terms.filt, gene %in% jgenes) %>%
  group_by(gene) %>%
  filter(rnk == min(rnk))
jgenes.terms <- make.names(jterms$termgene)

jsub <- as.data.frame(imput.mat[, jterms$term])
colnames(jsub) <- jgenes.terms

dat.exprs.proj <- data.frame(cell = rownames(imput.mat), jsub, stringsAsFactors = FALSE)
dat.merge.proj.exprs <- left_join(dat.merge.proj, dat.exprs.proj)


# pdf("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs/correct_var_LDA/k4me1_LDA_downstream_var_corrected.pdf")

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

# redo louvain
m.louv.after.redo <- ggplot(dat.merge.proj, aes(x = umap1, y = umap2, color = louvain_AfterCorrection)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

print(m.louv.after)

for (jgene in jgenes.terms){
  print(jgene)
  m <- PlotXYWithColor(dat.merge.proj.exprs, xvar = "umap1", yvar = "umap2", cname = jgene)
  print(m)
}
# dev.off()

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

# outf <- file.path("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/rdata", "k4me1_LDA_downstream_var_correted.RData")
# save(outf, topics.mat, terms.mat, rot.mat, topics.mat.proj, terms.mat.proj, annots.out, dat.merge, dat.merge.proj, file = outf)




# Load bulk ---------------------------------------------------------------

# Load bulk ---------------------------------------------------------------

inf.blood <- "/home/jyeung/hpc/scChiC/public_data/E-MTAB-3079-query-results.fpkms.tsv"

dat <- fread(inf.blood, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  group_by(Gene_Name, CellType) %>%
  summarise(FPKM = sum(FPKM)) %>%
  rowwise() %>%
  mutate(logFPKM = log2(FPKM + 1))

# normalize across samples?
ggplot(dat.long, aes(x = CellType, y = logFPKM)) + geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dat.mat <- tidyr::spread(dat.long %>%
                           ungroup() %>%
                           # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>%
                           mutate(gene = Gene_Name) %>%
                           dplyr::select(gene, CellType, logFPKM),
                         key = CellType, value = logFPKM)  %>%
  as.data.frame()
rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL

cnames.tmp <- colnames(dat.mat)
rnames.tmp <- rownames(dat.mat)
dat.mat <- preprocessCore::normalize.quantiles(as.matrix(dat.mat), copy = TRUE)  # strong normalization,
colnames(dat.mat) <- cnames.tmp
rownames(dat.mat) <- rnames.tmp

boxplot(dat.mat)

dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>%
  as.data.frame()
rownames(dat.norm.zscore.mat) <- dat.norm.zscore.mat$gene
dat.norm.zscore.mat$gene <- NULL


# Compare mystery islands to blood types  ---------------------------------


head(dat.merge.proj)

# find blood types frdom PCA

dat.umap.long.withPC <- left_join(dat.umap.long, data.frame(cell = rownames(topics.mat.proj), topics.mat.proj))

ggplot(dat.umap.long.withPC, aes(x = umap1, y = umap2, color = PC2)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.long.withPC, aes(x = umap1, y = umap2, color = PC4)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Find gene loadings of PC2

terms.dat <- data.frame(term = rownames(terms.mat.proj), PC2 = terms.mat.proj[, 1], PC3 = terms.mat.proj[, 2], PC4 = terms.mat.proj[, 3], stringsAsFactors = FALSE)

terms.dat <- left_join(terms.dat, annots.out$terms.filt) %>%
  group_by(term) %>%
  filter(weight == max(weight)) %>% 
  arrange(desc(PC4))


set.seed(0)
ngenes <- 50
genes.vec <- unique(terms.dat$gene)
immune.genes <- unique(terms.dat$gene[1:ngenes])
# immune.genes <- sample(genes.vec, size = ngenes, replace = FALSE)  # random

print(immune.genes)

jsub <-  subset(dat.norm.long, gene %in% immune.genes)

ggplot(jsub, aes(x = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = mean, .desc = TRUE), y = zscore)) + 
  geom_boxplot() + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# fit linear model to check if it's significant

jfit.sat <- lm(formula = zscore ~ celltype, data = jsub)
jfit.null <- lm(formula = zscore ~ 1, data = jsub)
jfit.anova <- anova(jfit.sat, jfit.null)
jfit.anova[["Pr(>F)"]][[2]]

# plot the whole thing as a density
ggplot(jsub, aes(x = zscore)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# permute zscores

ggplot(jsub, aes(x = exprs)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


Npermutes <- 1000
exprs.vec.lst <- lapply(seq(Npermutes), FUN = function(N){
  random.genes <- sample(genes.vec, size = ngenes, replace = FALSE)
  exprs.vec <- subset(dat.norm.long, gene %in% random.genes, select = c(celltype, gene, exprs))
  exprs.vec$iter <- N
  return(exprs.vec)
})

exprs.vec.lst.df <- exprs.vec.lst %>%
  bind_rows() %>%
  ungroup() %>%
  mutate(gene.select = "random",
         iter = as.character(iter))

# add my random one 
exprs.vec.lst.df <- bind_rows(exprs.vec.lst.df, subset(jsub, select = c(gene, celltype, exprs)) %>% mutate(iter = "FromPCA", gene.select = "FromPCA")) %>% 
  ungroup()

ggplot(exprs.vec.lst.df, aes(x = exprs, group = iter, fill = gene.select)) + geom_density()

# check medians?
exprs.vec.lst.df.sum <- exprs.vec.lst.df %>%
  # filter(gene.select == "random") %>%
  group_by(iter, gene.select) %>%
  summarise(exprs = median(exprs))

vline <- subset(exprs.vec.lst.df.sum, gene.select == "FromPCA")$exprs

ggplot(exprs.vec.lst.df.sum %>% filter(gene.select == "random"), aes(x = exprs, group = gene.select, fill = gene.select)) + geom_density() + geom_vline(xintercept = vline)

# ggplot(subset(dat.norm.long, gene %in% random.genes), aes(x = exprs)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



