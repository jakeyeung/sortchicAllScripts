# Jake Yeung
# Date of Creation: 2020-01-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/3-mix_BM_with_intestine.R
# Can we find immune signature by mixing? 




rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(irlba)



# Constants ---------------------------------------------------------------

binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * 1000

jsystem <- "Intestines"

jmark <- "k4me1"
jmark.bm <- "H3K4me1"

# jmark <- "k4me3"
# jmark.bm <- "H3K4me3"

jcutoff.ncuts.var <- 0.05


# Set up  -----------------------------------------------------------------

inf <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2020-01-11.10000_5000/lda_outputs.mat.Scraped.AllMerged.", jmark, ".TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_50000_25000.2020-01-11.K-30.binarize.FALSE/ldaOut.mat.Scraped.AllMerged.", jmark, ".TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_50000_25000.2020-01-11.K-30.Robj")
load(inf, v=T)
out.lda.int <- out.lda
count.mat.int <- count.mat

inf.bm <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.", jmark.bm, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark.bm, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
assertthat::assert_that(file.exists(inf.bm))
load(inf.bm, v=T)
out.lda.bm <- out.lda
count.mat.bm <- count.mat


# find mystery cells in intestine

topics.mat <- posterior(out.lda.int)$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  geom_hline(yintercept = -4.5)




dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "-"))
print(unique(dat.umap.long$experi))


immune.cells <- subset(dat.umap.long, umap2 < -4.5)$cell
lgr5.cells <- subset(dat.umap.long, grepl("stemcells", cell))$cell
dat.umap.long$gut.immune <- sapply(dat.umap.long$cell, function(x) x %in% immune.cells)
dat.umap.long$gut.lgr5 <- sapply(dat.umap.long$cell, function(x) grepl("stemcells", x))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = gut.immune)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  geom_hline(yintercept = -4.5)

bins.common <- intersect(rownames(count.mat.int), rownames(count.mat.bm))

count.mat.merge <- cbind(count.mat.int[bins.common, ], count.mat.bm[bins.common, ])

lsi.out <- scchicFuncs::RunLSI(as.matrix(count.mat.merge))

umap.out <- umap(lsi.out$u, config = jsettings)
dat.umap.long.lsi <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)

# dat.umap.long.lsi <- DoLouvain(lsi.out$u, jsettings, dat.umap.long.lsi)

dat.umap.long.lsi <- dat.umap.long.lsi %>%
  rowwise() %>% 
  mutate(experi = ifelse(startsWith(cell, "HVG"), "HVG", "PZ"),
         gut.immune = cell %in% immune.cells)

dat.umap.long.lsi <- dat.umap.long.lsi %>%
  rowwise() %>%
  mutate(gut.lgr5 = cell %in% lgr5.cells)


ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~experi)

ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)  + 
  ggtitle(paste(jmark, "BoneMarrow + Intestines PCA together shows mystery island is immune cells"))

ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2, color = gut.immune, shape = experi)) + geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)  + 
  ggtitle(paste(jmark, "BoneMarrow + Intestines PCA together shows mystery island is immune cells"))



# Show variance  ----------------------------------------------------------

tm.result.int <- posterior(out.lda.int)
tm.result.bm <- posterior(out.lda.bm)


dat.impute.log.int <- log2(t(tm.result.int$topics %*% tm.result.int$terms))
dat.impute.log.bm <- log2(t(tm.result.bm$topics %*% tm.result.bm$terms))

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var.int <- CalculateVarAll(dat.impute.log.int, jchromos)
dat.var.bm <- CalculateVarAll(dat.impute.log.bm, jchromos)

dat.var.merge <- rbind(dat.var.int, dat.var.bm)

dat.umap.long.lsi.var <- left_join(dat.umap.long.lsi, dat.var.merge)

ggplot(dat.umap.long.lsi.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~experi) + scale_color_viridis_c(direction = -1) + 
  ggtitle(paste(jmark, "BoneMarrow + Intestines overlap at regions of low variance"))

ggplot(dat.umap.long.lsi.var, aes(x = umap1, y = umap2, color = gut.lgr5)) + geom_point(alpha = 0.25) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  ggtitle(paste(jmark, "Lgr5 cells merge with bone marrow cells"))








