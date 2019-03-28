# Jake Yeung
# Date of Creation: 2019-03-23
# File: ~/projects/scchic/scripts/scripts_analysis/make_tables/compare_louvains.R
# Rerun LDA and compare louvain clustering


rm(list=ls())

library(ggplot2)
library(ggrepel)

library(dplyr)
library(hash)

library(umap)
library(igraph)
library(topicmodels)
library(JFuncs)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/GetMetaCellHash.R")


JaccardIndex <- function(x, y){
  numer <- length(intersect(x, y))
  denom <- length(union(x, y))
  return(numer / denom)
}


outdir <- "~/data/scchic/tables/bamlist_for_merging"
dir.create(outdir)

load("~/data/scchic/robjs/TFactivity_genelevels_objects.RData", v=T)

jmarks.all <- c("H3K4me1"="H3K4me1", "H3K4me3"="H3K4me3", "H3K27me3"="H3K27me3", "H3K9me3"="H3K9me3")

# out.lda <- ChooseBestLDA(out.lda)
# tm.result <- posterior(out.lda)
# topics.mat <- tm.result$topics

# jmark <- "H3K4me1"
out.lda.new.lst <- list()
inf1 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-15_20_25_30_35.Robj")
load(inf1, v=T)
out.lda.new.lst[["H3K4me1"]] <- ChooseBestLDA(out.lda)
inf2 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-15_20_25_30_35.Robj")
load(inf2, v=T)
out.lda.new.lst[["H3K4me3"]] <- ChooseBestLDA(out.lda)
inf3 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj")
load(inf3, v=T)
out.lda.new.lst[["H3K27me3"]] <- ChooseBestLDA(out.lda)
inf4 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
load(inf4, v=T)
out.lda.new.lst[["H3K9me3"]] <- ChooseBestLDA(out.lda)

topics.mat.new.lst <- lapply(out.lda.new.lst, function(out.lda){
  return(posterior(out.lda)$topics)
})


barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))


# topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
# topics.mat.lst <- lapply(tm.result.lst, function(tm.result) tm.result$topics)
dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks.all
topics.mat.lst <- lapply(tm.result.lst, function(tm.result) tm.result$topics)




# Do louvains for clustering  ------------------------------------------------------------

nn.louv <- c(27, 27, 60, 60)


jmetric.louv='euclidean' 
jmindist.louv=0.3
jseed.louv=123
custom.settings.louv.lst <- lapply(nn.louv, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))


clstr.hash.lst <- mapply(function(topics.mat, custom.settings.louv) DoLouvain(topics.mat, custom.settings.louv, dat.umap.long = NULL), topics.mat.lst, custom.settings.louv.lst, SIMPLIFY = FALSE)
# assign cluster to dat umap
dat.umap.long.lst <- lapply(jmarks.all, function(jmark){
  dat.umap <- dat.umap.lst[[jmark]]
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr.hash.lst[[jmark]][[x]]))
  dat.umap.long$mark <- jmark
  return(dat.umap.long)
})

  
  # nn.louv.new <- c(28, 35, 33, 27)
  nn.louv.new <- c(28, 35, 33, 31)
  jmindist.new <- c(0.2, 0.15, 0.3, 0.3)
  nn.new <- c(40, 30, 45, 27)
  # custom.settings.new.lst <- lapply(nn.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))
  custom.settings.new.lst <- mapply(function(x, y) GetUmapSettings(x, jmetric.louv, y, jseed.louv), nn.new, jmindist.new, SIMPLIFY = FALSE)
  custom.settings.louv.new.lst <- lapply(nn.louv.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))
  
  
  # custom.settings.new <- GetUmapSettings(nn=nn.louv[[1]], jmetric=jmetric.louv, jmindist=jmindist.louv, seed=jseed.louv)
  
  # do UMAP on new settings
  dat.umap.new.lst <- mapply(function(topics.mat, custom.settings) umap(topics.mat, config = custom.settings), 
                             topics.mat.new.lst, custom.settings.new.lst,
                             SIMPLIFY = FALSE)
  # do Louvain on new.louv settings
  clstr.hash.new.lst <- mapply(function(topics.mat, custom.settings) DoLouvain(topics.mat, custom.settings, dat.umap.long = NULL), 
                               topics.mat.new.lst, custom.settings.louv.new.lst)
  
  #  assign cluster to umap
  dat.umap.long.new.lst <- lapply(jmarks.all, function(jmark){
    dat.umap <- dat.umap.new.lst[[jmark]]
    dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
    dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr.hash.new.lst[[jmark]][[x]]))
    print(head(dat.umap.long))
    dat.umap.long$mark <- jmark
    return(dat.umap.long)
  })
  
  
  i <- 4
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#654321")
  ggplot(dat.umap.long.new.lst[[jmarks.all[[i]]]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) + scale_color_manual(values = cbPalette) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(nn.new[[i]], nn.louv.new[[i]]))
  

jmark <- "H3K4me3"
jmark <- "H3K4me1"
jmark <- "H3K27me3"
jmark <- "H3K9me3"

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#654321")
m1 <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=2.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark) + scale_color_manual(values = cbPalette)
# print(m1)
m1.new <- ggplot(dat.umap.long.new.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=2.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark) + scale_color_manual(values = cbPalette)
# print(m1.new)
multiplot(m1, m1.new, cols = 2)

# Compare cell clustering before and after

# save objects for future loading
save(dat.umap.long.new.lst, topics.mat.new.lst, out.lda.new.lst, barcodes, experihash, cellhash.bc, custom.settings.new.lst, custom.settings.louv.new.lst, file = "~/data/scchic/robjs/gene_levels_build95.Rdata")

jsize <- 1.2
jalpha <- 1
pdf("~/Dropbox/scCHiC_figs/FIG4_BM/primetime_plots/umaps_louvain_2019-03-24_550_counts_per_cell.pdf", useDingbats = FALSE)

jmark <- "H3K4me1"
m1 <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1)
m1.new <- ggplot(dat.umap.long.new.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1.new)
multiplot(m1, m1.new, cols = 2)
(frac.eryth <- nrow(subset(dat.umap.long.new.lst[[jmark]], louvain %in% c(1, 4))) / nrow(dat.umap.long.new.lst[[jmark]]))
dat.before <- dat.umap.long.lst[[jmark]]
dat.after <- dat.umap.long.new.lst[[jmark]]
# filter out common cells
cells.common <- intersect(dat.before$cell, dat.after$cell)
dat.before <- subset(dat.before, cell %in% cells.common)
dat.after <- subset(dat.after, cell %in% cells.common)
eryth.before <- subset(dat.before, louvain == 6)$cell
eryth.after <- subset(dat.after, louvain %in% c(4, 1))$cell
(eryth.ji <- JaccardIndex(eryth.before, eryth.after))

neutro.before <- subset(dat.before, louvain %in% c(3, 9))$cell
neutro.after <- subset(dat.after, louvain %in% c(8, 7))$cell
(neutro.ji <- JaccardIndex(neutro.before, neutro.after))

bcell.before <- subset(dat.before, louvain %in% c(1, 2, 8))$cell
bcell.after <- subset(dat.after, louvain %in% c(2, 3, 5))$cell
(bcell.ji <- JaccardIndex(bcell.before, bcell.after))
# 
# # no progenitors
# bcell.before <- subset(dat.before, louvain %in% c(1, 2))$cell
# bcell.after <- subset(dat.after, louvain %in% c(3, 5))$cell
# (bcell.ji <- JaccardIndex(bcell.before, bcell.after))

# bcell progenitors
# bcell.before.prog <- subset(dat.before, louvain %in% c(8))$cell
# bcell.after.prog <- subset(dat.after, louvain %in% c(2))$cell
# (bcell.ji.prog <- JaccardIndex(bcell.before, bcell.after))

tcell.before <- subset(dat.before, louvain %in% c(5))$cell
tcell.after <- subset(dat.after, louvain %in% c(9))$cell
(tcell.ji <- JaccardIndex(tcell.before, tcell.after))

nk.before <- subset(dat.before, louvain %in% c(4))$cell
nk.after <- subset(dat.after, louvain %in% c(10))$cell
(nk.ji <- JaccardIndex(nk.before, nk.after))

# fract erythryos



jmark <- "H3K4me3"
m1 <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1)
m1.new <- ggplot(dat.umap.long.new.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1.new)
multiplot(m1, m1.new, cols = 2)

(frac.eryth <- nrow(subset(dat.umap.long.new.lst[[jmark]], louvain %in% c(4))) / nrow(dat.umap.long.new.lst[[jmark]]))

# compare cells in the two clusters
dat.before <- dat.umap.long.lst[[jmark]]
dat.after <- dat.umap.long.new.lst[[jmark]]
# filter out common cells
cells.common <- intersect(dat.before$cell, dat.after$cell)
dat.before <- subset(dat.before, cell %in% cells.common)
dat.after <- subset(dat.after, cell %in% cells.common)
eryth.before <- subset(dat.before, louvain == 6)$cell
eryth.after <- subset(dat.after, louvain %in% c(4))$cell
(eryth.ji <- JaccardIndex(eryth.before, eryth.after))

neutro.before <- subset(dat.before, louvain %in% c(1, 3, 4))$cell
neutro.after <- subset(dat.after, louvain %in% c(1, 3, 6))$cell
(neutro.ji <- JaccardIndex(neutro.before, neutro.after))

bcell.before <- subset(dat.before, louvain %in% c(2, 5))$cell
bcell.after <- subset(dat.after, louvain %in% c(2, 5))$cell
(bcell.ji <- JaccardIndex(bcell.before, bcell.after))



jmark <- "H3K27me3"
m1 <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(paste(jmark, "after"), "after")) + scale_color_manual(values = cbPalette)
# print(m1)
m1.new <- ggplot(dat.umap.long.new.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1.new)
multiplot(m1, m1.new, cols = 2)
(frac.eryth <- nrow(subset(dat.umap.long.new.lst[[jmark]], louvain == 1)) /  nrow(subset(dat.umap.long.new.lst[[jmark]])))

# compare cells in the two clusters
dat.before <- dat.umap.long.lst[[jmark]]
dat.after <- dat.umap.long.new.lst[[jmark]]
# filter out common cells
cells.common <- intersect(dat.before$cell, dat.after$cell)
dat.before <- subset(dat.before, cell %in% cells.common)
dat.after <- subset(dat.after, cell %in% cells.common)


eryth.before <- subset(dat.before, louvain == 1)$cell
eryth.after <- subset(dat.after, louvain %in% c(1))$cell
(eryth.ji <- JaccardIndex(eryth.before, eryth.after))

neutro.before <- subset(dat.before, louvain %in% c(5, 6))$cell
neutro.after <- subset(dat.after, louvain %in% c(5, 6))$cell
(neutro.ji <- JaccardIndex(neutro.before, neutro.after))

bcell.before <- subset(dat.before, louvain %in% c(4, 7))$cell
bcell.after <- subset(dat.after, louvain %in% c(2, 7))$cell
(bcell.ji <- JaccardIndex(bcell.before, bcell.after))

# middle cells 
midcells.before <- subset(dat.before, louvain %in% c(3))$cell
midcells.after <- subset(dat.after, louvain %in% c(3))$cell
(midcells.ji <- JaccardIndex(midcells.before, midcells.after))


jmark <- "H3K9me3"
m1 <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1)
m1.new <- ggplot(dat.umap.long.new.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size=jsize, alpha = jalpha) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, "after")) + scale_color_manual(values = cbPalette)
# print(m1.new)
multiplot(m1, m1.new, cols = 2)
(frac.eryth <- nrow(subset(dat.umap.long.new.lst[[jmark]], louvain %in% c(4))) / nrow(dat.umap.long.new.lst[[jmark]]))

# compare cells in the two clusters
dat.before <- dat.umap.long.lst[[jmark]]
dat.after <- dat.umap.long.new.lst[[jmark]]
# filter out common cells
cells.common <- intersect(dat.before$cell, dat.after$cell)
dat.before <- subset(dat.before, cell %in% cells.common)
dat.after <- subset(dat.after, cell %in% cells.common)


eryth.before <- subset(dat.before, louvain == 6)$cell
eryth.after <- subset(dat.after, louvain %in% c(4))$cell
(eryth.ji <- JaccardIndex(eryth.before, eryth.after))

neutro.before <- subset(dat.before, louvain %in% c(2, 1, 7))$cell
neutro.after <- subset(dat.after, louvain %in% c(3, 5))$cell
(neutro.ji <- JaccardIndex(neutro.before, neutro.after))

bcell.before <- subset(dat.before, louvain %in% c(3, 4))$cell
bcell.after <- subset(dat.after, louvain %in% c(1, 2))$cell
(bcell.ji <- JaccardIndex(bcell.before, bcell.after))

# progens? cells
midcells.before <- subset(dat.before, louvain %in% c(5))$cell
midcells.after <- subset(dat.after, louvain %in% c(6))$cell
(midcells.ji <- JaccardIndex(midcells.before, midcells.after))

dev.off()


# Write to file  ----------------------------------------------------------

dat.umap.long.new.merged <- bind_rows(dat.umap.long.new.lst)

data.table::fwrite(dat.umap.long.new.merged, file = "~/data/scchic/tables/Cell_to_Louvain_550_counts_per_cell_2019-03-24.txt", sep = "\t", col.names = TRUE)


# Check no empty wells ----------------------------------------------------

indx.all <- seq(384)
hascell.indx <- c(seq(1:356),seq(360:379)+360)
empty.indx <- setdiff(indx.all, hascell.indx)

cell.indx <- sort(unique(as.numeric(sapply(dat.umap.long.new.merged$cell, function(x) strsplit(x, "cell")[[1]][[2]]))))


setdiff(empty.indx, cell.indx)


