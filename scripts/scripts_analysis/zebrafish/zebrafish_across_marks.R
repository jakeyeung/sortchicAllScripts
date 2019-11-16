# Jake Yeung
# Date of Creation: 2019-11-15
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_across_marks.R
# Across marks, nonenriched only?

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(JFuncs)

library(DESeq2)
library(preprocessCore)
library(CCA)

# Functions ---------------------------------------------------------------

DoUmapAndLouvain <- function(topics.mat, jsettings){
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
  dat.umap.long <- DoLouvain(topics.mat = topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long)
  return(dat.umap.long)
}

CollapseMatByLouvains <- function(count.mat, dat.umap.longs){
  jmat <- left_join(melt(as.matrix(count.mat)), dat.umap.longs %>% dplyr::select(c(cell, louvmark)), by = c("Var2" = "cell")) %>%
    group_by(louvmark, Var1) %>%
    summarise(count = sum(value)) %>%
    dcast(data = ., Var1 ~ louvmark) %>%
    as.data.frame()
  rownames(jmat) <- jmat$Var1
  jmat$Var1 <- NULL
  return(jmat)
}



# Some constants ----------------------------------------------------------

winsize <- 100000L
winsize <- 50000L
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
Kstrs <- c("30", "30", "30_35_50", "30")

infs <- mapply(function(jmark, Kstr) paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-ZFWKM-", 
                                    jmark, ".winsize_", winsize, ".merged.K-", 
                                    Kstr, ".binarize.FALSE/ldaOut.PZ-ChIC-ZFWKM-", 
                                    jmark, ".winsize_", winsize, ".merged.K-", Kstr, ".Robj"), jmarks, Kstrs)
assertthat::assert_that(all(file.exists(infs)))
names(infs) <- jmarks
names(jmarks) <- jmarks

# Load LDA outputs to get clustering --------------------------------------

lda.lst <- lapply(infs, function(inf){
  load(inf, v=T)
  if (length(out.lda) > 1){
    return(list(out.lda = out.lda[[1]], count.mat = count.mat))
  } else {
    return(list(out.lda = out.lda, count.mat = count.mat))
  }
})

# show umaps

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

topics.mats <- lapply(lda.lst, function(x) posterior(x$out.lda)$topics)
count.mats <- lapply(lda.lst, function(x) x$count.mat)
# filter common rows
common.rows <- Reduce(intersect, lapply(count.mats, function(x) rownames(x)))
count.mats <- lapply(count.mats, function(x) x[common.rows, ])


# Get pseudobulks from louvain  -------------------------------------------

dat.umap.longs <- lapply(topics.mats, DoUmapAndLouvain, jsettings)

dat.umap.longs <- lapply(jmarks, function(jmark){
  return(dat.umap.longs[[jmark]] %>% mutate(mark = jmark))
}) %>% bind_rows() %>%
  mutate(louvain = paste("louvain", louvain, sep = "_"),
         louvmark = paste(louvain, mark, sep = "-"))
  
m.umaps <- ggplot(dat.umap.longs, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~mark) + ggtitle(paste("Winsize:", winsize)) 

count.mats.collapsed.lst <- lapply(count.mats, CollapseMatByLouvains, dat.umap.longs)
count.mats.collapsed <- do.call(cbind, count.mats.collapsed.lst)

jmarks.sub <- c("H3K4me1", "H3K27me3")
jgrp <- paste(jmarks.sub, collapse = "|")

# inmat <- log2(t(count.mats.collapsed[, grepl(jgrp, colnames(count.mats.collapsed))]))
# jcutoff <- 6
# genes.keep <- which(apply(inmat, 2, min) > jcutoff) 
# 
# dat.pca.long <- data.frame(louvmark = rownames(dat.pca$x), pc1 = dat.pca$x[, 1], pc2 = dat.pca$x[, 2], stringsAsFactors = FALSE) %>%
#   rowwise() %>%
#   mutate(louvain = strsplit(strsplit(louvmark, "-")[[1]][[1]], "\\.")[[1]][[2]],
#          mark = strsplit(louvmark, "-")[[1]][[2]])
# 
# cbPalette <- c("#ff9f7d", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#eb9d01", "#7fbedf", "#696969")
# m.umaps.filt <- ggplot(subset(dat.umap.longs, mark %in% jmarks.sub), aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#   facet_wrap(~mark) + ggtitle(paste("Winsize:", winsize)) + 
#   scale_color_manual(values = cbPalette)
# m.pca <- ggplot(dat.pca.long, aes(x = pc1, y = pc2, shape = mark, color = louvain, label = louvmark)) + geom_point(size = 5) + geom_text_repel() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
#   scale_color_manual(values = cbPalette)
# 
# multiplot(m.pca, m.umaps.filt, cols = 2)


# Get correlations from raw data ------------------------------------------

# correlate H3K4me1 versus H3K4me3 first, easy??

source("/Users/yeung/projects/scchic/scripts/Rfunctions/IntegrateData.R")

X <- as.matrix(count.mats.collapsed.lst[[jmarks.sub[[1]]]])
Y <- as.matrix(count.mats.collapsed.lst[[jmarks.sub[[2]]]])

# normalize count size
ctype.metadata <- data.frame(ctype = c(colnames(X), colnames(Y)), mark = c(rep(jmarks.sub[[1]], ncol(X)), rep(jmarks.sub[[2]], ncol(Y))))
rownames(ctype.metadata) <- ctype.metadata$ctype
dds <- DESeqDataSetFromMatrix(countData = cbind(X, Y), colData = ctype.metadata, design = ~ mark)
vsd <- vst(dds)

XY.norm <- normalize.quantiles(assay(vsd))
rownames(XY.norm) <- rownames(assay(vsd)); colnames(XY.norm) <- colnames(assay(vsd))
# scale?
# XY.norm <- t(scale(t(XY.norm), center = TRUE, scale = TRUE))

X.norm <- XY.norm[, grepl(jmarks.sub[[1]], colnames(XY.norm))]
Y.norm <- XY.norm[, grepl(jmarks.sub[[2]], colnames(XY.norm))]

# # scale
# X.norm <- t(scale(t(X.norm), center = TRUE, scale = TRUE))
# Y.norm <- t(scale(t(Y.norm), center = TRUE, scale = TRUE))

X.norm.vst <- assay(vsd)[, grepl(jmarks.sub[[1]], colnames(assay(vsd)))]
Y.norm.vst <- assay(vsd)[, grepl(jmarks.sub[[2]], colnames(assay(vsd)))]


# Do PCA ------------------------------------------------------------------

dat.pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
dat.pca <- prcomp(t(XY.norm), center = TRUE, scale. = TRUE)

plot(dat.pca)

plot(dat.pca$x[, 1], dat.pca$x[, 2])
text(dat.pca$x[, 1], dat.pca$x[, 2], labels = rownames(dat.pca$x))

dat.pca.long <- data.frame(louvmark = rownames(dat.pca$x), pc1 = dat.pca$x[, 1], pc2 = dat.pca$x[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(louvain = strsplit(louvmark, "-")[[1]][[1]], 
         mark = strsplit(louvmark, "-")[[1]][[2]])

cbPalette <- c("#ff9f7d", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#eb9d01", "#7fbedf", "#696969")
m.umaps.filt <- ggplot(subset(dat.umap.longs, mark %in% jmarks.sub), aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  facet_wrap(~mark) + ggtitle(paste("Winsize:", winsize)) + 
  scale_color_manual(values = cbPalette)
m.pca <- ggplot(dat.pca.long, aes(x = pc1, y = pc2, shape = mark, color = louvain, label = louvmark)) + geom_point(size = 5) + geom_text_repel() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette)

multiplot(m.pca, m.umaps.filt, cols = 2)


# Do CCA ------------------------------------------------------------------

cca.out <- CCA::cc(X.norm, Y.norm)
# cca.out <- CCA::cc(X, Y)
# cca.out <- cancor(X.norm, Y.norm)

# # cca.out <- jCanonCor(mat1 = X, mat2 = Y, l2.norm = TRUE, k = 2, use.irlba = FALSE)
# cca.out <- jCanonCor(mat1 = X.norm, mat2 = Y.norm, l2.norm = FALSE, k = 2, use.irlba = FALSE)
# # cca.out <- jCanonCor(mat1 = X.norm.vst, mat2 = Y.norm.vst, l2.norm = TRUE, k = 2, use.irlba = FALSE)
# 
# mat3 <- t(scale(X.norm, center = TRUE, scale = TRUE)) %*% scale(Y.norm, center = TRUE, scale = TRUE)
# mat3 <- matcor(X.norm, Y.norm)
# cca.svd <- svd(mat3$XYcor)
# cca.svd$u <- cca.svd$u / sqrt(sum(cca.svd$u ^ 2))
# cca.svd$v <- cca.svd$v / sqrt(sum(cca.svd$v ^ 2))
# cca.out <- list(u = cca.svd$u, v = cca.svd$v, d = cca.svd$d)

cca.data <- rbind(cca.out$xcoef, cca.out$ycoef)
colnames(x = cca.data) <- paste0("CC", 1:ncol(cca.data))

# icca.data) <- c(colnames(X.norm), colnames(Y.norm))
# rownames(cca.data) <- c(colnames(X.norm), colnames(Y.norm))
# 
# # flip signs for aesthetics?
# cca.data.flip <- apply(cca.data, MARGIN = 2, function(x){
#   if(sign(x[1]) == -1) {
#     x <- x * -1
#   }
#   return(x)
# })


# get loadings by projecting
# how are these feature loadings calculated? Just projection
# embeds <- cca.data.flip
embeds <- cca.data
# average across embeds
# loads <- t(t(embeds) %*% rbind(t(X.norm), t(Y.norm)))
# loads <- t(t(embeds) %*% rbind(t(X.norm), t(Y.norm)))
# take min function?
embeds1 <- cca.out$u
embeds2 <- cca.out$v

cca1 <- 1
cca2 <- 2
plot(embeds[, cca1], embeds[, cca2], pch = 20)
text(embeds[, cca1], embeds[, cca2],
     labels = rownames(embeds))

dat.cca.long <- data.frame(louvmark = rownames(embeds), cca1 = embeds[, 1], cca2 = embeds[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(louvain = strsplit(louvmark, "-")[[1]][[1]], 
         mark = strsplit(louvmark, "-")[[1]][[2]])
  
m.cca <- ggplot(dat.cca.long, aes(x = cca1, y = cca2, label = louvmark, shape = mark, color = louvain)) +
  geom_point(size = 5) + geom_text_repel() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

multiplot(m.cca, m.umaps.filt, cols = 2)

# plot correlations 

# ggpairs(as.data.frame(cbind(X.norm, Y.norm)))


# infer using marker genes  -----------------------------------------------


# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
inf.WKM <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"
dat.bulk <- readRDS(inf.WKM)
dat.bulk.mat <- dcast(subset(dat.bulk, select = c(gene, celltype, exprs)), gene ~ celltype, value.var = "exprs")
rownames(dat.bulk.mat) <- dat.bulk.mat$gene; dat.bulk.mat$gene <- NULL

# set up reference data
zscore.cutoff <- "top_100"
if (is.character(zscore.cutoff)){
  rnk.keep <- as.numeric(strsplit(zscore.cutoff, "_")[[1]][[2]])
  assertthat::assert_that(is.numeric(rnk.keep))
  dat.bulk.keep <- dat.bulk %>%
    group_by(celltype) %>%
    mutate(rnk = rank(-zscore)) %>%
    filter(rnk < rnk.keep)
} else {
  dat.bulk.keep <- dat.bulk %>%
    group_by(gene) %>%
    filter(max(abs(zscore)) > zscore.cutoff)
}

rnames.common <- Reduce(intersect, lapply(count.mats, function(x) rownames(x)))
genes.chic <- sapply(rnames.common, function(x) strsplit(x, ";")[[1]][[2]])
ref.genes.keep <- intersect(as.character(dat.bulk.keep$gene), genes.chic)

print(length(ref.genes.keep))

# For each set of genes, plot distirbution of read counts across l --------

# eryth genes
# jctype <- "erythrocytes"
print(unique(dat.bulk.keep$celltype))

jctype <- "erythrocytes"
jctype <- "monocytes"
jctype <- "HSPCs"
jctype <- "lymphocytes"
jctype <- "neutrophils"
jctype <- "thrombocytes"

jgenes <- as.character(subset(dat.bulk.keep, celltype == jctype)$gene)

print(head(sort(jgenes)))

# take h3k4me1 across louvains
X.norm.filt <- CollapseRowsByGene(X.norm, as.long = TRUE) %>%
  dplyr::rename(gene = Var1, louvain = Var2) %>%
  group_by(gene) %>%
  mutate(zscore = scale(count, center = TRUE, scale = TRUE))
m.ctypes <- ggplot(subset(X.norm.filt, gene %in% jgenes), aes(x = louvain, y = zscore)) + geom_boxplot() + geom_jitter(width = 0.25, height = 0)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jctype)
multiplot(m.ctypes, m.umaps.filt, cols = 2)






