# Jake Yeung
# Date of Creation: 2020-01-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/find_bad_cells.R
# Find bad cells


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)
library(DropletUtils)

# Load raw ----------------------------------------------------------------

inf <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2020-01-11.10000_5000/lda_outputs.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_50000_25000.2020-01-11.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_50000_25000.2020-01-11.K-30.Robj"
load(inf, v=T)


tm.result <- posterior(out.lda)
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

# Select high deviance genes  ---------------------------------------------

ncuts.dat <- data.frame(cell = colnames(count.mat), ncuts = colSums(count.mat))
nvec <- ncuts.dat$ncuts
pvec <- rowSums(count.mat) / sum(count.mat)

bd.gene <- apply(X = count.mat, MARGIN = 1, FUN = function(xvec){
  binomial_deviance(xvec, sum(xvec) / sum(nvec), nvec)
})

plot(sort(bd.gene, decreasing = TRUE))


# Select high deviance cells  ---------------------------------------------

ncuts.dat.gene <- data.frame(gene = rownames(count.mat), ncuts.gene = rowSums(count.mat))
nvec.gene <- ncuts.dat.gene$ncuts
bd.cell <- apply(X = count.mat, MARGIN = 2, FUN = function(xvec){
  binomial_deviance(xvec, sum(xvec) / sum(nvec.gene), nvec.gene)
})


# multinom
pseudocount <- 1
prop <- min(colSums(count.mat)) / colSums(count.mat)
count.mat.downsamp <- DropletUtils::downsampleMatrix(count.mat, prop)

md <- apply(X = count.mat, MARGIN = 2, FUN = function(jcol){
  # jcol <- jcol + pseudocount
  pnull <- 1 / length(jcol)
  multinomial_deviance(jcol, pnull)
})

# coltest <- rep(1, ncol(count.mat))
# coltest[1] <- 2
# coltest[2] <- 0
# pnull <- rep(1 / length(coltest), length(coltest))
# multinomial_deviance(coltest, pnull)


md.dat <- data.frame(cell = names(md), mdev = md, stringAsFactor = FALSE)
plot(md)

# look at bad cell in raw data


jcell <- as.character(subset(dat.var, cell.var.within.sum.norm == max(cell.var.within.sum.norm))$cell)
jcell <- as.character(subset(dat.var, cell.var.within.sum.norm == min(cell.var.within.sum.norm))$cell)
subset(dat.merge, cell == jcell)

jcells <- colnames(count.mat)
names(jcells) <- jcells
names(jchromos) <- jchromos
rows.i.lst <- lapply(jchromos, function(jchromo){
  terms.keep <- grep(jchromo, rownames(count.mat))
})

mat.sum.by.chromo <- lapply(rows.i.lst, function(terms.keep){
  return(colMeans(count.mat[terms.keep, ]))
  # return(colSums(count.mat[terms.keep, ]))
})

mat.sum.by.chromo.merged <- do.call(rbind, mat.sum.by.chromo)

reads.by.chromo <- lapply(jcells, function(jcell){
  reads.by.chromo <- mat.sum.by.chromo.merged[, jcell]
  reads.by.chromo.dat <- data.frame(chromo = names(reads.by.chromo), ncuts = unlist(reads.by.chromo), stringsAsFactors = FALSE) %>%
    mutate(cell = jcell)
  return(reads.by.chromo.dat)
}) %>% 
  bind_rows()

# get inter-chromosomal variance?
reads.by.chromo.sum <- reads.by.chromo %>% 
  group_by(cell) %>%
  filter(!chromo %in% c("chrX", "chrY")) %>%
  summarise(intrachr.var = var(ncuts)) %>%
  left_join(., dat.var) %>%
  left_join(., ncuts.dat) %>%
  mutate(log2var = log2(intrachr.var))

reads.by.chromo.sum.frompca <- reads.by.chromo %>% 
  group_by(cell) %>%
  filter(chromo %in% c("chr11", "chr19", "chrX", "chrY")) %>%
  mutate(high.chr = ifelse(chromo %in% c("chr11", "chr19"), TRUE, FALSE)) %>%
  group_by(cell, high.chr) %>%
  summarise(ncuts = mean(ncuts)) %>%
  group_by(cell) %>% 
  summarise(intrachr.var = var(ncuts)) %>%
  left_join(., dat.var) %>%
  left_join(., ncuts.dat) %>%
  mutate(log2var = log2(intrachr.var))

ggplot(reads.by.chromo.sum.frompca, aes(x = intrachr.var, y = cell.var.within.sum.norm)) + geom_point()
ggplot(reads.by.chromo.sum.frompca, aes(x = ncuts, y = intrachr.var)) + geom_point()


ggplot(reads.by.chromo.sum, aes(x = intrachr.var, y = cell.var.within.sum.norm)) + geom_point()
ggplot(reads.by.chromo.sum, aes(x = log2(intrachr.var), y = cell.var.within.sum.norm, color = ncuts)) + geom_point()+ scale_color_viridis_c()
ggplot(reads.by.chromo.sum, aes(x = cell.var.across, y = intrachr.var)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(reads.by.chromo.sum, aes(x = ncuts, y = intrachr.var, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1)
ggplot(reads.by.chromo.sum, aes(x = ncuts, y = log2(intrachr.var), color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1)

# do loess fit?
fit.var <- loess(log2var ~ ncuts, reads.by.chromo.sum)

fit.dat <- data.frame(cell = reads.by.chromo.sum$cell, log2var = reads.by.chromo.sum$log2var, ncuts = reads.by.chromo.sum$ncuts, log2var.pred = predict(fit.var, newdata = reads.by.chromo.sum)) %>%
  left_join(., dat.var)

ggplot(fit.dat, aes(x = log2var - log2var.pred, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(fit.dat, aes(x = log2var - log2var.pred, y = ncuts, color = cell.var.within.sum.norm)) + geom_point() + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)




# plot test cells 
jcell <- as.character(subset(dat.var, cell.var.within.sum.norm == min(cell.var.within.sum.norm))$cell)
jcell <- as.character(subset(dat.var, cell.var.within.sum.norm == max(cell.var.within.sum.norm))$cell)

jsub <- subset(reads.by.chromo.sum, ncuts > quantile(ncuts, 0.99)) %>% arrange(desc(intrachr.var))
# jsub <- subset(reads.by.chromo.sum, ncuts > quantile(ncuts, 0.99)) %>% arrange(intrachr.var)
jcell <- subset(jsub, intrachr.var == min(intrachr.var))$cell[[1]]
jcell <- subset(jsub, intrachr.var == max(intrachr.var))$cell[[1]]

ggplot(subset(reads.by.chromo, cell == jcell), aes(x = chromo, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jcell)


# PCA on counst mat
# X <- data.table::dcast(reads.by.chromo %>% filter(!chromo %in% c("chrX", "chrY")), chromo ~ cell, value.var = "ncuts")
X <- data.table::dcast(reads.by.chromo, chromo ~ cell, value.var = "ncuts")
X.mat <- t(as.matrix(subset(X, select = -chromo)))
colnames(X.mat) <- X$chromo
X.mat <- scale(X.mat, center = TRUE, scale = FALSE)
X.pca <- prcomp(X.mat, center = FALSE, scale. = FALSE)

dat.pca <- data.frame(cell = rownames(X.pca$x), pc1 = X.pca$x[, 1], pc2 = X.pca$x[, 2], stringsAsFactors = FALSE) %>%
  left_join(., dat.var) %>%
  left_join(., reads.by.chromo.sum)

ggplot(dat.pca, aes(x = pc1, y = pc2, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1)
ggplot(dat.pca, aes(x = pc1, y = pc2, color = intrachr.var)) + geom_point() + scale_color_viridis_c(direction = -1)
# ggplot(X.pca, aes(x = X.pca$x[, 1], X.pca$x[, 2]) + geom_point()

# what is pc2?
dat.pca.loadings <- data.frame(chromo = rownames(X.pca$rotation), X.pca$rotation)

ggplot(dat.pca.loadings, aes(x = PC1, y = PC2, label = chromo)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_text()

# show a cell with ihg PC2 loading
jcell <- subset(dat.pca, pc2 == max(pc2))$cell

ggplot(subset(reads.by.chromo, cell == jcell), aes(x = chromo, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jcell)


jsd <- signif(sd(unlist(reads.by.chromo)), digits = 2)
plot(unlist(reads.by.chromo), main = paste(jcell, "\n", jsd))


plot(sort(bd.cell, decreasing = TRUE))
bd.cell.dat <- data.frame(cell = names(bd.cell), bdev = bd.cell, stringsAsFactors = FALSE)


dat.merge <- left_join(bd.cell.dat, dat.var)
dat.merge <- left_join(dat.merge, ncuts.dat)
dat.merge  <- left_join(dat.merge, md.dat)

ggplot(dat.merge, aes(x = bdev, y = cell.var.within.sum.norm, color = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = mdev, y = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge, aes(x = mdev, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = bdev, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggplot(dat.merge, aes(x = mdev, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + scale_y_log10()
ggplot(dat.merge, aes(x = ncuts, y = log(bdev))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge, aes(x = ncuts, y = log(bdev), color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# bdev is function of ncuts, regress it out?
dat.merge.sort <- dat.merge %>% arrange(ncuts)
bdev.loess.out <- loess(formula = bdev ~ ncuts, data = dat.merge.sort)
# mdev.loess.out <- loess(formula = bdev ~ ncuts, data = dat.merge.sort)

plot(dat.merge.sort$ncuts, dat.merge.sort$bdev)
points(dat.merge.sort$ncuts, predict(loess.out), type = "l", col = 'blue', lwd = 10)

dat.merge.sort$bdev.pred <- predict(loess.out, newdata = dat.merge.sort)

dat.merge.sort <- dat.merge.sort %>%
  rowwise() %>%
  mutate(bdev.norm = bdev - bdev.pred)

ggplot(dat.merge.sort, aes(x = bdev.norm, y = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggplot(dat.merge.sort, aes(x = bdev.norm, y = cell.var.within.sum.norm, color = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c()
ggplot(dat.merge.sort, aes(x = bdev.norm, y = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.merge.sort, aes(x = ncuts, y = bdev.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1)

dat.merge.sort <- dat.merge.sort %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))


ggplot(dat.merge.sort, aes(x = bdev.norm, y = cell.var.within.sum.norm, color = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  facet_wrap(~plate)

ggplot(dat.merge.sort, aes(y = bdev.norm, color = cell.var.within.sum.norm, x = ncuts)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~plate)
  

ggplot(dat.merge, aes(x = bdev / ncuts, y = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.merge, aes(x = ncuts, y = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)




