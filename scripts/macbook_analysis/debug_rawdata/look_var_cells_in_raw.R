# Jake Yeung
# Date of Creation: 2020-01-17
# File: ~/projects/scchic/scripts/macbook_analysis/debug_rawdata/look_var_cells_in_raw.R
# Look at low variance cells in the raw data

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(ggrepel)


# Functions ---------------------------------------------------------------




# Show variance -----------------------------------------------------------

jmark <- "H3K4me3"
# inf <- "/Users/yeung/data/scchic/from_cluster/LDA_outputs/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/LDA_outputs/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
load(inf, v=T)

tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"),
         batch = ClipLast(cell, jsep = "-"))

# Show plate effects?? ----------------------------------------------------

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~experi)

# check bad topic, is it Dock2? 
# annot.out <- AnnotateBins2(terms.mat, top.thres = 0.995, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")


# # bind bad topic? 
# cells.linneg <- grep("Linneg", rownames(topics.mat), value = TRUE)
# jsub <- topics.mat[cells.linneg, ]
# 
# colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")
# 
# plot(colSums(jsub))

# dat.umap.topics <- left_join(dat.umap.long, data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE))
# 
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_30")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_21")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_14")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_26")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_26")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_22")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_6")
# PlotXYWithColor(dat.umap.topics, xvar = "umap1", yvar = "umap2", cname = "topic_19")
# 
# terms.sub <- annot.out$terms.annot
# 
# subset(terms.sub, grepl("Dock2", termgene))
# subset(terms.sub, grepl("Zfp982", termgene))
# subset(terms.sub, grepl("Pak3", termgene))
# subset(terms.sub, grepl("Grik2", termgene))
# subset(terms.sub, grepl("Hlf", termgene))

# can we see this in the raw data ?  --------------------------------------

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.meta <- data.frame(cell = colnames(count.mat),
                       ncuts = colSums(count.mat),
                       sparsity = apply(count.mat, 2, function(x) nnzero(x) / length(x)), 
                       stringsAsFactors = FALSE) %>%
  left_join(., dat.var)

ggplot(dat.meta, aes(x = ncuts, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.meta, aes(x = sparsity, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.meta, aes(x = cell.var.across, y = cell.var.within.sum.norm, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c(direction = -1)

# cell var across should be robust even in single cells??

count.devres <- null_residuals(m = as.matrix(count.mat), mod = "binomial")


dat.var.raw <- CalculateVarAll(count.devres, jchromos)

count.mat.norm <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")

# sum across chromos?
# jcell <- "PZ-ChIC-Bl6-BM-H3K4me1-Index2-12-09-19_355"
# xvec <- count.mat.norm[, jcell]
# 
# names(jchromos) <- jchromos
# jsum <- lapply(jchromos, function(jchromo){
#   terms.keep <- grepl(paste0("^", jchromo, ":"), names(xvec))
#   jsub <- xvec[terms.keep]
#   return(sum(jsub))
# })

jsum <- SumAcrossChromos(count.mat.norm, jchromos, colfunction = mean)

jsum.var <- jsum %>%
  filter(!chromo %in% c("chrX", "chrY")) %>%
  group_by(cell) %>%
  summarise(ncuts.var.chromo = var(log2(ncuts))) %>%
  # summarise(ncuts.var.chromo = var(ncuts)) %>%
  left_join(., dat.meta, by = "cell") %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"),
         batch = ClipLast(cell, jsep = "-"))

ggplot(jsum.var, aes(x = ncuts.var.chromo, y = cell.var.within.sum.norm)) + geom_point()
ggplot(jsum.var, aes(x = , y = cell.var.within.sum.norm)) + geom_point()

ggplot(jsum.var, aes(x = ncuts.var.chromo, y = cell.var.across, color = cell.var.within.sum.norm)) + geom_point() + 
  scale_x_log10() + scale_y_log10()  + scale_color_viridis_c(direction = -1)

ggplot(jsum.var, aes(x = ncuts.var.chromo, y = cell.var.within.sum.norm, color = log10(ncuts))) + geom_point() +
  scale_color_viridis_c()

ggplot(jsum.var, aes(x = ncuts.var.chromo, color = cell.var.within.sum.norm, y = ncuts)) + geom_point() +
  scale_color_viridis_c(direction = -1)  + scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsum.var, aes(x = ncuts, color = cell.var.within.sum.norm, y = cell.var.within.sum.norm)) + geom_point() +
  scale_color_viridis_c(direction = -1)  + scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

ggplot(jsum.var, aes(x = ncuts.var.chromo, color = cell.var.within.sum.norm, y = log10(ncuts))) + geom_point() +
  scale_color_viridis_c(direction = -1)

ggplot(jsum.var, aes(x = ncuts.var.chromo, color = cell.var.within.sum.norm, y = log10(ncuts))) + geom_point() +
  scale_color_viridis_c(direction = -1)

ggplot(jsum.var, aes(x = ncuts.var.chromo)) + geom_density() + 
  scale_color_viridis_c(direction = -1)

ggplot(jsum.var, aes(x = ncuts.var.chromo, y = ncuts)) + geom_point() + 
  scale_color_viridis_c(direction = -1) + scale_x_log10() + scale_y_log10()

# related to sparsity probably





