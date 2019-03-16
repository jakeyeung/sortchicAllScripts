# Jake Yeung
# Date of Creation: 2019-02-14
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/plot_marker_genes_all_4_marks_more_bins.R
# 

rm(list=ls())

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(umap)
library(ggrepel)
library(biomaRt)
library(igraph)  # louvain
library(Gviz)
library(GenomicRanges)


source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load data ---------------------------------------------------------------

inmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"

meanfilt <- 10
# jbin <- "FALSE"; kstr <- "15_20_25_30"
jbin <- "TRUE"; kstr <- "15_20_25_30_35"

# jmark <- "H3K27me3"
jmark <- "H3K9me3"
# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
indir <- paste0("lda_outputs.meanfilt_", 
                meanfilt,
                ".cellmin_100.cellmax_500000.binarize.", jbin, 
                ".no_filt")
fname <- paste0("lda_out_meanfilt.BM-", jmark, ".CountThres0.K-", kstr, ".Robj")

inf <- file.path(inmain, indir, fname)

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

out.lda <- ChooseBestLDA(out.lda)
(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

# settings for UMAP
nn=30
nnterms <- 15
jmetric='euclidean' 
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], -dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))




top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

subset(top.peaks, topic == 2)

# subset(top.peaks, topic == 26)
# subset(top.peaks, topic == 27)
# subset(top.peaks, topic == 21)
# subset(top.peaks, topic == 13)

regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd), 
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)

top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))

print(subset(top.peaks.annotated, grepl("Hox", SYMBOL)))
print(subset(top.peaks.annotated, grepl("Sox6", SYMBOL)))
print(subset(top.peaks.annotated, grepl("Hbb", SYMBOL)))
print(subset(top.peaks.annotated, grepl("Hba", SYMBOL)))

subset(top.peaks.annotated, topic == 2)

# translate beta to log fold change?
mat.norm <- t(tm.result$topics %*% tm.result$terms)  # this should give normalized signal, without the poisson noise?

mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")

# visualize a peak
jpeak <- "chr12:115560000-115660000"

# get clstring based on topic2 vs all
# jlab <- as.numeric(dat.umap$layout[, 1] < -3)

dat.umap.louv <- dat.umap

cell.indx <- hash(rownames(dat.umap.louv$knn$indexes), dat.umap.louv$knn$indexes[, 1])
cell.indx.rev <- hash(dat.umap.louv$knn$indexes[, 1], rownames(dat.umap.louv$knn$indexes))
nr <- nrow(dat.umap.louv$knn$indexes)
nc <- ncol(dat.umap.louv$knn$indexes)
edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
colnames(edgelist) <- c("from", "to")
for (vertex.i in seq(nr)){
  istart <- nc*(vertex.i - 1)+1
  iend <- nc*vertex.i
  edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
  edgelist[istart : iend, 2] <- sapply(dat.umap.louv$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
  # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
}
g <- graph_from_data_frame(edgelist, directed=FALSE)
g.out <- cluster_louvain(g, weights = NULL)
V(g)$color <- g.out$membership
clstr <- hash(g.out$names, g.out$membership)

# clstr <- hash(rownames(dat.umap$layout), jlab)

jchromo <- "chr12"
# get data for chromo 7
jpeaks <- grep(jchromo, rownames(mat.norm), value = TRUE)
x <- as.data.frame(mat.norm[jpeaks, ])
x.long <- data.frame(exprs = unlist(x), cell = rep(colnames(x), each = nrow(x)),
                     coord = rep(rownames(x), ncol(x)), stringsAsFactors = FALSE)
x.long$louvain.orig <- sapply(x.long$cell, function(x) clstr[[x]])
x.long$louvain <- x.long$louvain.orig
x.long$exprs <- x.long$exprs * 10^6


dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
# dat.umap.long <- left_join(dat.umap.long, topics.mat)

dat.umap.long$louvain <- sapply(as.character(dat.umap.long$cell), function(x) clstr[[x]])

# m.louvain <- ggplot(dat.umap.long %>% mutate(louvain = ifelse(louvain == "4", TRUE, FALSE)), aes(x = umap1, y = umap2, color = as.character(louvain))) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.louvain)

jstart <- as.numeric(GetStart(jpeak)) -  5 * 10 ^ 5
jend <- as.numeric(GetEnd(jpeak)) + 5 * 10 ^ 5
PlotGTrack(x.long %>% mutate(louvain = ifelse(louvain == "4", TRUE, FALSE)),
           jstart, 
           jend, 
           mart.obj, gen = "mm10", chr = jchromo, jheight = 1.5)


# Load all 4 datasets -----------------------------------------------------


