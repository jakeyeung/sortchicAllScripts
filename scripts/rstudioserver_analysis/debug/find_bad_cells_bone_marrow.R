# Jake Yeung
# Date of Creation: 2020-01-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/find_bad_cells.R
# Find bad cells

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)
library(DropletUtils)


inf <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
load(inf, v=T)

tm.result <- posterior(out.lda)
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
names(jchromos) <- jchromos
rows.i.lst <- lapply(jchromos, function(jchromo){
  terms.keep <- grep(jchromo, rownames(count.mat))
})

mat.sum.by.chromo <- lapply(rows.i.lst, function(terms.keep){
  return(colMeans(count.mat[terms.keep, ]))
  # return(colSums(count.mat[terms.keep, ]))
})

mat.sum.by.chromo.merged <- do.call(rbind, mat.sum.by.chromo)

jcells <- colnames(count.mat)
names(jcells) <- jcells
reads.by.chromo <- lapply(jcells, function(jcell){
  reads.by.chromo <- mat.sum.by.chromo.merged[, jcell]
  reads.by.chromo.dat <- data.frame(chromo = names(reads.by.chromo), ncuts = unlist(reads.by.chromo), stringsAsFactors = FALSE) %>%
    mutate(cell = jcell)
  return(reads.by.chromo.dat)
}) %>% 
  bind_rows()

# is there signal ? 

X <- data.table::dcast(reads.by.chromo, chromo ~ cell, value.var = "ncuts")
X.mat <- t(as.matrix(subset(X, select = -chromo)))
colnames(X.mat) <- X$chromo
X.mat <- scale(X.mat, center = FALSE, scale = FALSE)
X.pca <- prcomp(X.mat, center = FALSE, scale. = FALSE)

dat.pca <- data.frame(cell = rownames(X.pca$x), pc1 = X.pca$x[, 1], pc2 = X.pca$x[, 2], stringsAsFactors = FALSE) %>%
  left_join(., dat.var) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"))
  # left_join(., reads.by.chromo.sum)

ggplot(dat.pca, aes(x = pc1, y = pc2, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1) 
ggplot(dat.pca, aes(x = pc1, y = pc2, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1)  + 
  facet_wrap(~experi)

# check loadings

# what is pc2?
dat.pca.loadings <- data.frame(chromo = rownames(X.pca$rotation), X.pca$rotation)

ggplot(dat.pca.loadings, aes(x = PC1, y = PC2, label = chromo)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_text()

# show a cell with ihg PC2 loading
jcell <- subset(dat.pca, pc2 == max(pc2))$cell
jcell <- subset(dat.pca, pc2 == min(pc2))$cell

jcell <- (dat.pca %>% arrange(pc2))$cell[[15]]
jcell <- (dat.pca %>% arrange(desc(pc2)))$cell[[1]]
jcell <- (dat.pca %>% arrange(pc2))$cell[[2]]

ggplot(subset(reads.by.chromo, cell == jcell), aes(x = chromo, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jcell)

# mean of chr11 and chr19 vs autosomes?
reads.by.chromo.sum.frompca <- subset(reads.by.chromo, !chromo %in% c("chrX", "chrY")) %>%
  mutate(high.chr = chromo %in% c("chr11", "chr19")) %>%
  group_by(cell, high.chr) %>%
  summarise(ncuts = mean(ncuts)) %>%
  group_by(cell) %>%
  summarise(ncuts.log2fc = ncuts[[2]] - ncuts[[1]]) %>%
  left_join(., dat.var)

ggplot(reads.by.chromo.sum.frompca, aes(x = ncuts.log2fc, y = cell.var.within.sum.norm)) + geom_point()


