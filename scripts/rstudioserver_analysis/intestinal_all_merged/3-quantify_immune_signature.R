# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/3-quantify_immune_signature.R
# See if we can qauntify immune signature in intestines 

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

# Load data ---------------------------------------------------------------


# inf.lda <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.2019-12-22.K-30.Robj"
inf.lda <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2020-01-11.10000_5000/lda_outputs.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_50000_25000.2020-01-11.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.k4me1.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_50000_25000.2020-01-11.K-30.Robj"

load(inf.lda, v=T)

tm.result <- posterior(out.lda)


topics.mat <- tm.result$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# plot sparsity? 
ncuts <- data.frame(cell = colnames(count.mat), ncuts = colSums(count.mat), stringsAsFactors = FALSE)
sparsity <- data.frame(cell = colnames(count.mat), frac.zeros = apply(count.mat, 2, function(x) 1 - nnzero(x) / length(x)), stringsAsFactors = FALSE)

dat.merge.withmeta <- left_join(dat.merge, ncuts)
dat.merge.withmeta <- left_join(dat.merge.withmeta, sparsity)

ggplot(dat.merge.withmeta, aes(x = umap1, y = umap2, color = frac.zeros)) + 
  geom_point(size = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = 1)

# ggplot(dat.merge.withmeta, aes(x = umap1, y = umap2, color = log10(ncuts))) + 
ggplot(dat.merge.withmeta, aes(x = umap1, y = umap2, color = ncuts)) + 
  geom_point(size = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = 1)

ggplot(dat.merge.withmeta, aes(x = ncuts, y = frac.zeros)) + geom_point() + scale_x_log10()


# Do naive fraction of counts and see if it's indicative of varian --------


count.mat.frac <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")

count.mat.frac[which(as.matrix(count.mat.frac) == 0)] <- NA

par(mfrow=c(2, 2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# plot a low var cell
(jsub <- subset(dat.merge.withmeta, cell.var.within.sum.norm == min(cell.var.within.sum.norm)))
jcell <- jsub$cell[[1]]
plot(density(count.mat.frac[, jcell], na.rm = TRUE), main = jsub$cell.var.within.sum.norm[[1]])
sort(unique(count.mat.frac[, jcell]))
sort(unique(count.mat.frac[, jcell]))
jtab <- table(count.mat.frac[, jcell])
plot(names(jtab), unlist(jtab), log = "y")
# fit a line?

jtab.dat <- data.frame(unlist(jtab)) %>%
  mutate(log10freq = log10(Freq),
         number = as.numeric(Var1))
(jfit <- lm(formula = log10freq ~ number, data = jtab.dat))

(jsub <- subset(dat.merge.withmeta, cell.var.within.sum.norm == max(cell.var.within.sum.norm)))
jcell <- jsub$cell[[1]]
plot(density(count.mat.frac[, jcell], na.rm = TRUE), main = jsub$cell.var.within.sum.norm[[1]])
sort(unique(count.mat.frac[, jcell]))
sort(unique(count.mat[, jcell]))
print(length(unique(sort(unique(count.mat[, jcell])))))
jtab <- table(count.mat.frac[, jcell])
plot(names(jtab), unlist(jtab), log = "y")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

jtab.dat <- data.frame(unlist(jtab)) %>%
  mutate(log10freq = log10(Freq),
         number = as.numeric(Var1))
(jfit <- lm(formula = log10freq ~ number, data = jtab.dat))


# sum across chromosomes ? 

names(jchromos) <- jchromos
dat.chromo.avg <- lapply(jchromos, function(jchromo){
  jrows <- grepl(paste0("^", jchromo), rownames(count.mat))
  chromo.avg <- data.frame(cell = colnames(count.mat), chromo.mean = colMeans(count.mat[jrows, ]), chromo = jchromo, stringsAsFactors = FALSE)
  return(chromo.avg)
}) %>%
  bind_rows()

dat.chromo.sum <- dat.chromo.avg %>%
  group_by(cell) %>%
  # mutate(cell.var.across.raw = stats::mad(log2(chromo.mean)))
  mutate(cell.var.across.raw = stats::mad(chromo.mean))

dat.chromo.sum <- left_join(dat.chromo.sum, dat.var)
dat.chromo.sum <- left_join(dat.chromo.sum, ncuts)

ggplot(dat.chromo.sum, aes(x = cell.var.within.sum.norm, y = cell.var.across.raw)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.chromo.sum, aes(x = cell.var.across.raw, y = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.chromo.sum, aes(x = ncuts, y = cell.var.across.raw)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10()
# ubset(d)

ggplot(dat.chromo.sum, aes(x = ncuts, y = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10()





