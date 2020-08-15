# Jake Yeung
# Date of Creation: 2020-08-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_glmpca_with_spikein.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

library(irlba)

library(hash)
library(igraph)
library(umap)

library(topicmodels)


make.plots <- TRUE

# jmark <- "H3K9me3"
# jmark <- "H3K4me1"
# jmark <- "H3K27me3"
jmark <- "H3K4me3"

# jsuffix <- ".penalty_5.by_plate.RData"
# jsuffix <- ".topn_5000.glmpcaout.penalty_5.by_plate.RData"
jsuffix <- ".glmpcaout.penalty_10RData"
jprefix <- ""

if (make.plots){
  pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562_cellcycle/dim_reduction_summaries_with_spikeins.", jmark, ".", Sys.Date(), jprefix, jsuffix, ".pdf")
  pdf(pdfout, useDingbats = FALSE)
}




# Load data  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_K562_spikein/lda_outputs.K562_count_tables_50000.", jmark, ".", jprefix, "K-30.binarize.FALSE/ldaOut.K562_count_tables_50000.", jmark, ".", jprefix, "K-30.Robj"))
assertthat::assert_that(file.exists(inf.lda))

# inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, ".G1_G2_S.glmpcaout", jsuffix))
inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, jprefix, jsuffix))
assertthat::assert_that(file.exists(inf.glmpca))


inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, ".", jprefix, "rds"))
inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.spike))


mat <- readRDS(inf)

load(inf.spike, v=T)

dat.spikeins.mat <- dat.spikeins.mat[colnames(mat), ]


totalcounts <- data.frame(cell = colnames(mat), totalcounts = colSums(mat), stringsAsFactors = FALSE)

spikeincounts.vec <- dat.spikeins.mat[colnames(mat), ]$spikeincounts
names(spikeincounts.vec) <- dat.spikeins.mat[colnames(mat), ]$samp

cellcycle.lst = list("0" = "0_G1", "1" = "1_S", "2" = "2_G2/M")
cellcycle.hash <- hash::hash(cellcycle.lst)

dat.spikeins.mat <- dat.spikeins.mat %>%
  rowwise() %>%
  mutate(cellid = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
         indx = strsplit(samp, "_")[[1]][[2]], 
         rowcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[1]],
         colcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[2]],
         is.empty = rowcoord <= 8 & colcoord == 1) %>%
  left_join(., totalcounts, by = c("samp" = "cell"))

dat.spikeins.mat$experi <- sapply(dat.spikeins.mat$samp, function(x) ClipLast(x, jsep = "_"))

lsi.out <- RunLSI2(count.mat = as.matrix(mat), tf.method = "spikeins", idf.method = "inverselog", n.components = 30, .log = FALSE, .center = FALSE, .truncated = TRUE, spikeincounts = spikeincounts.vec)

dim(lsi.out$u)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap.lsi <- DoUmapAndLouvain(lsi.out$u, jsettings) %>%
  left_join(., dat.spikeins.mat, by = c("cell" = "samp"))


ggplot(dat.umap.lsi, aes(x = umap1, y = umap2, color = experi)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~experi) + 
  ggtitle("LSI")
  
ggplot(dat.umap.lsi, aes(x = umap1, y = umap2, color = log2(totalcounts / spikeincounts))) + geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LSI")

ggplot(dat.umap.lsi, aes(x = umap1, y = umap2, color = log2(totalcounts / spikeincounts))) + geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LSI") + 
  facet_wrap(~experi)

ggplot(dat.umap.lsi, aes(x = experi, y = log2(totalcounts / spikeincounts))) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LSI")

# Load LDA ----------------------------------------------------------------

load(inf.lda, v=T)

dat.umap.lda <- DoUmapAndLouvain(posterior(out.lda)$topics, jsettings = jsettings) %>%
  left_join(., dat.spikeins.mat, by = c("cell" = "samp"))

ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi) + 
  ggtitle("LDA")

ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LDA")

ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = log2(totalcounts / spikeincounts))) + 
  geom_point() +
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LDA")

ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = log10(totalcounts))) + 
  geom_point() +
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LDA")

# Load GLMPCA -------------------------------------------------------------

load(inf.glmpca, v=T)

dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings = jsettings) %>%
  left_join(., dat.spikeins.mat, by = c("cell" = "samp"))

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi) + 
  ggtitle("GLMPCA")

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = log2(totalcounts / spikeincounts))) + 
  geom_point() +
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("GLMPCA")

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = log2(totalcounts / spikeincounts))) + 
  geom_point() +
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("GLMPCA")



# Add variance ?  ---------------------------------------------------------

tm.result <- posterior(out.lda)

dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))

jchromos <- paste("chr", c(seq(19)), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)


ggplot(dat.umap.lda %>% left_join(., dat.var), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() +
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LDA")

ggplot(dat.umap.lsi %>% left_join(., dat.var), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() +
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("LSI")

ggplot(dat.umap.glmpca %>% left_join(., dat.var), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() +
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("GLMPCA")

if (make.plots){
  dev.off()
  
}
