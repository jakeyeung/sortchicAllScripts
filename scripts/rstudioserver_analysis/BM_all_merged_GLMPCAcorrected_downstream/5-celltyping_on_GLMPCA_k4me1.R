# Jake Yeung
# Date of Creation: 2020-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/5-celltyping_on_GLMPCA.R
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

library(scchicFuncs)


# Functions ---------------------------------------------------------------


# Load data  --------------------------------------------------------------


jmark <- "H3K4me1"
print(jmark)

jexperi <- "AllMerged"

mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))
load(inf.glm, v=T)


inf.glmpca.annot <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
assertthat::assert_that(file.exists(inf.glmpca.annot))
load(inf.glmpca.annot, v=T)

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

inf.annot <- "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
load(inf.annot, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))


# Load annotations  -------------------------------------------------------

# findloading for Neutrophils for examle
unique(dat.umap.glm.fillNAs$cluster)


jcluster <- "Eryth_topic27"
jcluster <- "Neutrophils_topic23"

jclust.vec <- as.numeric(dat.umap.glm.fillNAs$cluster == jcluster)
jclust.vec[is.na(jclust.vec)] <- 0

assertthat::assert_that(all(rownames(glm.out$factors) == dat.umap.glm.fillNAs$cell))
glm.out$factors

jcors <- apply(glm.out$factors, MARGIN = 2, FUN = function(jcol){
  cor(jclust.vec, jcol, method = "spearman")
})

dat.umap.glm.merged <- left_join(dat.umap.glm.fillNAs, data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE)) %>%
  rowwise() %>%
  mutate(cond = GetCondFromSamp(cell, mark = jmark))

(jdim <- names(jcors)[which.max(abs(jcors))])
print(jcors)
barplot(jcors)

# plot top 
ggplot(dat.umap.glm.merged, aes_string(x = "umap1", y = "umap2", color = jdim)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = 1) + facet_wrap(~cond)

sum(apply(glm.out$factors[, c(1,2)], MARGIN = 1, prod))

sum(apply(glm.out$loadings[, c(1,2)], MARGIN = 1, prod))

# SVD an example
pca.out <- prcomp(glm.out$factors, center = TRUE, scale. = FALSE)
sum(apply(pca.out$rotation[, c(1, 2)], MARGIN = 1, prod))


# Compare loadings with public data ---------------------------------------

# get common rows then plot correlations?

inf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
annot.out <- AnnotateCoordsFromList(coords.vec = rownames(glm.out$loadings), inf.tss = inf.tss)

# x <- glm.out$loadings[, jdim]
# rename region_coord to symbol?

rnames.old <- rownames(glm.out$loadings)
coord2gene <- hash::hash(annot.out$regions.annotated$region_coord, annot.out$regions.annotated$SYMBOL)

glm.out$loadings[, jdim]

loadings.renamed <- data.frame(coord = rownames(glm.out$loadings), 
                               gene = sapply(rownames(glm.out$loadings), AssignHash, jhash = coord2gene),
                               glm.out$loadings, stringsAsFactors = FALSE)

# add gene expression (zscore) to the loadings?
loadings.sub <- loadings.renamed[, c("coord", "gene", jdim)]
# add zscore across cells
dat.sum.zscore <- reshape2::dcast(data = dat.sum.long, formula = gene ~ celltype, value.var = "zscore")
# loadings.sub <- left_join(loadings.sub, dat.sum.zscore)

#  get common genesD
genes.common <- intersect(dat.sum.zscore$gene, loadings.sub$gene)

dat.sum.zscore.merge <- inner_join(dat.sum.zscore, loadings.sub) %>%
  filter(!duplicated(gene))

ctypes <- colnames(dat.sum.norm)

for (ctype in ctypes){
  print(ctype)
  x <- dat.sum.zscore.merge[[jdim]]
  y <- dat.sum.zscore.merge[[ctype]]
  jcor <- cor(x, y, method = "spearman")
  plot(x, y, xlab = jdim, ylab = ctype, pch = 20, main = paste(ctype, jcor))
}

# take top N genes and plot? 
jtop <- 150
jsub.sorted <- loadings.sub[order(loadings.sub[[jdim]], decreasing = FALSE), ]
genes.keep <- jsub.sorted$gene[1:jtop]

# plot boxplots of filtered genes
library(forcats)

ggplot(dat.sum.long %>% filter(gene %in% genes.keep) %>% mutate(celltype = fct_reorder(.f = celltype, .x = zscore, .desc=TRUE)), aes(x = celltype, y = zscore)) + geom_boxplot()   


# Plot imputed expression of marker genes? --------------------------------

dat.impute <- t(as.matrix(glm.out$factors) %*% as.matrix(t(glm.out$loadings)))

# jgene <- "S100a7"
jgene <- "Bach2"
jgene <- "Siglech"
jgene <- "Irf4"

jgene <- "Hlf"

jgene <- "Ltf"
(jsub <- subset(annot.out$regions.annotated, grepl(jgene, SYMBOL)))
jterm <- jsub$region_coord[[1]]

exprs.dat <- data.frame(cell = colnames(dat.impute), exprs = dat.impute[jterm, ], stringsAsFactors = FALSE)

dat.umap.glm.impute <- left_join(dat.umap.glm.fillNAs, exprs.dat)

# if we average across many regions maybe we can get something interesting? 
ggplot(dat.umap.glm.impute, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
  scale_color_viridis_c(direction = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(jgene, jterm))


# Use topic loadings from K27me3 to plot GLMPCA outputs --------------------

# load active mark

jmark.act <- "H3K4me1"
jmark.repress <- "H3K27me3"
inf.act <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark.act, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark.act, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
assertthat::assert_that(file.exists(inf.act))

inf.repress <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark.repress, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark.repress, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
assertthat::assert_that(file.exists(inf.repress))

load(inf.act, v=T)
count.mat.act <- count.mat
out.lda.act <- out.lda

tm.result.act <- posterior(out.lda)
rownames(tm.result.act$terms) <- paste("topic", rownames(tm.result.act$terms), sep = "")
colnames(tm.result.act$topics) <- paste("topic", colnames(tm.result.act$topics), sep = "")


load(inf.repress, v=T)
count.mat.repress <- count.mat
out.lda.repress <- out.lda

tm.result.repress <- posterior(out.lda)
rownames(tm.result.repress$terms) <- paste("topic", rownames(tm.result.repress$terms), sep = "")
colnames(tm.result.repress$topics) <- paste("topic", colnames(tm.result.repress$topics), sep = "")


# Get topics from K27me3 and ask what the data looks like on K4me1? -------


jtop.repress <- "topic30"

jtop.repress <- "topic1"
jtop.repress <- "topic16"


jtop.repress <- "topic28"
jtop.repress <- "topic30"

jtop.repress <- "topic8"

jtop.repress <- "topic17"
jtop.repress <- "topic22"
jtop.repress <- "topic26"
jtop.repress <- "topic29"
jtop.repress <- "topic9"

jtop.repress <- "topic30"

jtop.repress <- "topic23"

jtop.repress <- "topic8"



jtop.repress <- "topic20"

jtop.repress <- "topic17"

jtop.repress <- "topic22"

jtop.repress <- "topic27"

jtop.repress <- "topic9"

jtop.repress <- "topic27"


jtop.repress <- "topic22"
# for (jtop.act in rownames(tm.result.act$terms)){
repress.loadings <- sort(tm.result.repress$terms[jtop.repress, ], decreasing = FALSE)
# act.loadings <- sort(tm.result.repress$terms["topic9", ], decreasing = TRUE)
(coords.keep.repress <- names(repress.loadings)[1:300])

coords.common <- intersect(coords.keep.repress, rownames(dat.impute))
print(length(coords.common))

# plot on K4me1 UMAP?
exprs.dat <- data.frame(cell = colnames(dat.impute), exprs = colMeans(dat.impute[coords.common, ]), stringsAsFactors = FALSE) %>%
  left_join(., dat.umap.glm.fillNAs)

ggplot(exprs.dat, mapping = aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(direction = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jtop.repress)


# Compare with GLMPCA -----------------------------------------------------

# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# dat.act.lda <- DoUmapAndLouvain(tm.result.act$topics, jsettings = jsettings)
# ggplot(dat.act.lda, aes(x = umap1, y = umap2)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 

# find loadings for neutrophil topic
jtop.act <- "topic14"  # HSCs
jtop.act <- "topic25"  # random
jtop.act <- "topic5"  # HSCs

jtop.act <- "topic25"  # HSCs


jtop.act <- "topic6"  # neutrophils
act.loadings <- sort(tm.result.act$terms[jtop.act, ], decreasing = TRUE)
(coords.keep.act <- names(act.loadings)[1:100])

rows.i <- which(rownames(dat.impute) %in% coords.keep.act)
# rows.i <- sample(x = seq(nrow(dat.impute)), size = 100, replace = FALSE)
print(rows.i)

# plot in GLMPCA, average across egions 
jsub <- dat.impute[rows.i, ]

# check genes assigned in jsub
jsub.annot <- subset(annot.out$regions.annotated, region_coord %in% rownames(jsub))

# print(jsub.annot$SYMBOL)

dat.sub <- data.frame(cell = colnames(dat.impute), exprs = colMeans(jsub), stringsAsFactors = FALSE) %>%
  left_join(., subset(dat.umap.glm.fillNAs, select = c(umap1, umap2, cell)))

print(head(dat.sub))

ggplot(dat.sub, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
  scale_color_viridis_c(direction = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jtop.act)



#colSums 
# 
# # plot(dat.sum.zscore.merge$dim1, dat.sum.zscore.merge$core, pch = 20, main = 'test')
# 
# 
# 
# 
# # plot zscore correlations???
# glm.out$loadings[1:5, 1:5]
# 
# 
# 
# 
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# umap.out <- umap(topics.mat, config = jsettings)
# dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# 
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)





