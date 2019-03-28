# Jake Yeung
# Date of Creation: 2019-03-22
# File: ~/projects/scchic/scripts/scripts_analysis/pseudobulk_analysis/louvain_pseudobulk_proportions_argument.R
# Similar proportions?


library(ggplot2)
library(ggrepel)

library(dplyr)
library(hash)

library(umap)
library(igraph)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

outdir <- "~/data/scchic/tables/bamlist_for_merging"
dir.create(outdir)

load("~/data/scchic/robjs/TFactivity_genelevels_objects.RData", v=T)

jmarks.all <- list("H3K4me1" = "H3K4me1", "H3K4me3" = "H3K4me3", "H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")


barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))




nn.louv <- c(27, 27, 60, 60)
jmetric.louv='euclidean' 
jmindist.louv=0.4
jseed.louv=123

custom.settings.louv.lst <- lapply(nn.louv, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))

# topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
topics.mat.lst <- lapply(tm.result.lst, function(tm.result) tm.result$topics)

dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks.all

topics.mat.lst <- lapply(tm.result.lst, function(tm.result) tm.result$topics)

clstr.hash.lst <- mapply(function(topics.mat, custom.settings.louv) DoLouvain(topics.mat, custom.settings.louv, dat.umap.long = NULL), topics.mat.lst, custom.settings.louv.lst, SIMPLIFY = FALSE)

# assign cluster to dat umap

dat.umap.long.lst <- lapply(jmarks.all, function(jmark){
  dat.umap <- dat.umap.lst[[jmark]]
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr.hash.lst[[jmark]][[x]]))
  dat.umap.long$mark <- jmark
  return(dat.umap.long)
})

# get bam to louvain list
jmark <- "H3K4me1"
jmark <- "H3K4me3"
jmark <- "H3K27me3"
jmark <- "H3K9me3"

print(head(dat.umap.long.lst[[jmark]]))


# Plot the 4 maps side by side  -------------------------------------------

dat.umap.long <- bind_rows(dat.umap.long.lst)

library(JFuncs)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
mlst <- lapply(dat.umap.long.lst, function(dat.umap.long){
  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 0.5) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values=cbPalette)
  return(m)
})
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 2)

# similar proportions?
eryth.clust <- list("H3K4me1" = 6, "H3K4me3" = 6, "H3K27me3" = 1, "H3K9me3" = 6)  # eryth
# eryth.clust <- list("H3K4me1" = 4, "H3K4me3" = 6, "H3K27me3" = 2, "H3K9me3" = 6)  # NK cells
# eryth.clust <- list("H3K4me1" = 3, "H3K4me3" = 5, "H3K27me3" = 5, "H3K9me3" = 2)  # NK cells
# eryth.clust <- list("H3K4me1" = c(1, 2, 4, 8, 5), "H3K4me3" = c(1, 3, 4), "H3K27me3" = c(4, 2, 7), "H3K9me3" = c(3, 4))  # NK cells
# eryth.clust <- list("H3K4me1" = c(3, 9, 7), "H3K4me3" = c(2, 5), "H3K27me3" = c(5, 6, 3), "H3K9me3" = c(7, 1, 2))  # NK cells
eryth.clust <- list("H3K4me1" = 4, "H3K4me3" = c(2, 5), "H3K27me3" = c(5, 6, 3), "H3K9me3" = c(7, 1, 2))  # NK cells

eryth.props <- mapply(function(mark, dat.umap.long){
  ntotal <- length(unique(dat.umap.long.lst[[mark]]$cell))
  jclst <- eryth.clust[[mark]]
  neryth <- length(unique(subset(dat.umap.long.lst[[mark]], louvain %in% jclst)$cell))
  return(neryth / ntotal)
}, jmarks.all, dat.umap.long.lst, SIMPLIFY = FALSE)
print(eryth.props)

dat.umap.long$indx <- sapply(dat.umap.long$cell, function(x) strsplit(x, "cell")[[1]][[2]])


# Look at empty wells -----------------------------------------------------

indx.all <- seq(384)
input.wells <- c(seq(1:356),seq(360:379)+360)
empty.wells <- setdiff(indx.all, input.wells)

dat.empty <- subset(dat.umap.long, indx %in% empty.wells)

empty.cells <- dat.empty$cell

mlst <- lapply(dat.umap.long.lst, function(dat.umap.long){
  dat.umap.long <- dat.umap.long %>% rowwise() %>% mutate(empty = cell %in% empty.cells)
  dat.umap.long <- RankOrder(dat.umap.long, cname = "empty", out.cname = "cell.lab.bank")
  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = empty)) + 
    geom_point(size = 2) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values=cbPalette)
  return(m)
})
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 2)

# what are the counts in empty wells?

names(count.mat.lst) <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

count.mat.lst.sums <- unlist(lapply(count.mat.lst, function(count.mat) colSums(count.mat)))


library(topicmodels)

library(dplyr)
library(ggplot2)
library(JFuncs)
library(topicmodels)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(hash)
library(umap)
library(scales)

jmarks <- jmarks.all
meanfilt <- 10
Kstr.nobin <- "15_20_25_30"
infs.nobin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)

cellsums.lst <- sort(unlist(lapply(jmarks, function(x) Matrix::colSums(out.objs.nobin[[x]]$count.mat))))
names(cellsums.lst) <- sapply(names(cellsums.lst), function(x) strsplit(x, "\\.")[[1]][[2]])

plot(density(cellsums.lst))

empty.vec <- sapply(names(cellsums.lst), function(x) ifelse(x %in% empty.cells, "Empty", ""))
col.vec <- sapply(empty.vec, function(x) ifelse(x == "", "blue", "red"))
cellsums.lst.filt <- cellsums.lst[which(empty.vec == "Empty")]

plot(x = seq(length(cellsums.lst)), y = cellsums.lst, col = "lightblue", pch = 20, cex = 5, log = "y")
points(x = seq(length(cellsums.lst.filt)), y = cellsums.lst.filt, col = "red", pch = 20, cex = 5)

# new cutoff
jcutoff <- 510

# how many cells do you lose if you do 510 gone?

dsub <- cellsums.lst[which(cellsums.lst > 510)]

dsub <- cellsums.lst[which(cellsums.lst > 550)]

# how many do you keep... PER MARK??

# how many are gone?
cells.gone <- setdiff(names(cellsums.lst), names(dsub))

cells.gone.by.mark <- sapply(cells.gone, function(x) strsplit(x, "_")[[1]][[2]])

table(cells.gone.by.mark)

mlst <- lapply(dat.umap.long.lst, function(dat.umap.long){
  dat.umap.long <- dat.umap.long %>% rowwise() %>% mutate(empty = cell %in% cells.gone)
  dat.umap.long <- RankOrder(dat.umap.long, cname = "empty", out.cname = "cell.lab.bank")
  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + 
    geom_point(size = 2) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values=cbPalette)
  return(m)
})
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 2)
# text(x = seq(length(cellsums.lst)), y = cellsums.lst, labels = empty.vec), cex = 2)

     