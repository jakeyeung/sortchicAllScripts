# Jake Yeung
# Date of Creation: 2019-01-24
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_LDA_25kb_H3K4me1.R
# Any difference with 25kb?

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

library(igraph)  # louvain

library(Gviz)
library(GenomicRanges)


source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

source("scripts/Rfunctions/PlotFunctions.R")

# Functions ---------------------------------------------------------------




# Load --------------------------------------------------------------------

jbin <- "TRUE"
jchip <- "H3K4me3"
jchip <- "H3K4me1"
# from t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/lda_analysis.h3k27me3.dropbox/lda_outputs.meanfilt_1000.merge_Windows2.cellmin_NA.cellmax_NA/lda_out.meanfilt.K-12.Robj
# from t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-5_10_15_20_25.Robj
# inf <- paste0("/private/tmp/lda_output_binned/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_10_15_20_25.Robj")
if (jbin){
  inf <- paste0("/tmp/ldaAnalysisBins_MetaCell_binsize_25kb/lda_outputs.binsize_25kb.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE/lda_out_meanfilt.BM-", jchip, "-25kb.CountThres0.K-15_20_25_30.Robj")
} else {
  stop("Not run yet")
  inf <- paste0("/private/tmp/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_15_25.Robj")
}

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# out.lda <- out.lda[[5]]
# print(out.lda@loglikelihood)
out.lda <- ChooseBestLDA(out.lda)
(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

nn=40
# nn=15
nnterms <- 15
jmetric='euclidean' 
jmindist=0.4
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# check your umap settings
# jpeak <- "chr7:103800000-103900000"
# PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)


# Plot dat umap -----------------------------------------------------------
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))

# Plot terms umap ---------------------------------------------------------
top.thres <- 0.995
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
top.regions <- unique(unlist(topic.regions))
terms.mat <- t(tm.result$terms)[top.regions, ]
# dat.umap.terms <- umap(terms.mat, config = custom.settings.terms)
# # downsample rows for plotting purposes
# downsamp.i <- sample(seq(nrow(dat.umap.terms$layout)), size = round(0.1 * nrow(dat.umap.terms$layout)), replace = FALSE)
# jcol.rgb.terms <- lapply(seq(kchoose), ColorsByGamma, terms.mat[downsamp.i, ], c("pink", "red", "darkred"))
# par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
# mapply(function(jcol.rgb.term, jtopic){
#   plot(dat.umap.terms$layout[downsamp.i, 1], dat.umap.terms$layout[downsamp.i, 2],
#        col = jcol.rgb.term, pch = 20, asp = 0.75,
#        main = paste("Peak Weights, T", jtopic))
# }, jcol.rgb.terms, seq(kchoose))

# how many cells in the island?
cell.assign <- apply(topics.mat, 1, which.max)

top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

# annotate regions?
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

hit.peaks <- subset(regions.annotated, abs(distanceToTSS) < 50000 & grepl("Hbb", SYMBOL))$region_coord

jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)

# how top hits for specific topics
# Progenitor cellsare in Topic 12?
topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)
dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)

# plot topic12 
jtopic <- 12
jcol.rgb <- jcol.rgbs[[jtopic]]
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# make it pretty
m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `12`)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
print(m)

# plot top hits for topic 12
print(subset(top.peaks.annotated, topic == 12), n = 50)
print(subset(top.peaks.annotated, topic == 10), n = 50)
print(subset(top.peaks.annotated, topic == 7), n = 50)

jgene <- c("Tgm2", "Mmrn1", "Trim48", "Pdzk1ip1", "Mllt3", "Mecom")
jgene <- c("Vldlr", "Uhrf1", "Rrm2", "Lig1", "Tipin")
jgene <- c("F2r", "Itga2b", "Zfp385a", "Zfpm1", "Plek", "Cd9", "Zeb2")

# translate beta to log fold change?
mat.norm <- t(tm.result$topics %*% tm.result$terms)
mat.norm <- mat.norm[top.regions, ]

# label using louvain clustering?


cell.indx <- hash(rownames(dat.umap$knn$indexes), dat.umap$knn$indexes[, 1])
cell.indx.rev <- hash(dat.umap$knn$indexes[, 1], rownames(dat.umap$knn$indexes))
nr <- nrow(dat.umap$knn$indexes)
nc <- ncol(dat.umap$knn$indexes)
# nc <- 4
# edgelist <- matrix(NA, nrow = nr * nc, ncol = 3)
edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
colnames(edgelist) <- c("from", "to", "weight")
for (vertex.i in seq(nr)){
  istart <- nc*(vertex.i - 1)+1
  iend <- nc*vertex.i
  edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
  edgelist[istart : iend, 2] <- sapply(dat.umap$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
  # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
}
g <- graph_from_data_frame(edgelist, directed=FALSE)
g.out <- cluster_louvain(g, weights = NULL)
V(g)$color <- g.out$membership
clstr <- hash(g.out$names, g.out$membership)

dat.umap.long$louvain <- sapply(dat.umap.long$cell, function(x) clstr[[x]])

# plot umap with louvain 
m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

# plot graph with edges?
# https://stackoverflow.com/questions/5364264/how-to-control-the-igraph-plot-layout-with-fixed-positions

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
coords <- layout.auto(g)
plot.igraph(simplify(g),
            layout = dat.umap$layout[V(g)$name, ],
            vertex.label = NA,
            edge.curved=FALSE,
            label = NA,
            edge.width = 0.5,
            vertex.size = 1)



# Merge cells and plot hits -----------------------------------------------

library(biomaRt)

gen <- "mm10"
chr <- "chr7"

jstart <- 115510106
jend <- 116999509
mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
itrack <- IdeogramTrack(genome=gen, chromosome=chr)
gtrack <- BiomartGeneRegionTrack(genome=gen, chromosome=chr, symbol = "Sox6", stacking = "full",
                       biomart = mart.obj,
                       name="Region",showId=T)
# data
print(subset(top.peaks.annotated, topic == 12), n = 50)
# out <- ParseCoord("chr7:115,468,872-116,040,744")
jpeak <- "chr7:115920000-116020000"
x <- mat.norm[jpeak, ]
x.long <- data.frame(exprs = x, cell = names(x), coord = jpeak, stringsAsFactors = FALSE)
x.long$louvain <- sapply(x.long$cell, function(x) clstr[[x]])
x.sum <- x.long %>%
  group_by(louvain, coord) %>%
  summarise(exprs = mean(exprs)) %>%
  mutate(start = as.numeric(GetStart(coord)),
         end = as.numeric(GetEnd(coord)),
         chromo = GetChromo(coord),
         data = exprs)
dat <- x.sum[which.max(x.sum$data), ]$exprs
jstarts <- rep(as.numeric(GetStart(jpeak)), 1)
jends <- rep(as.numeric(GetEnd(jpeak)), 1)
dtrack <- DataTrack(data=dat, 
                    start = jstarts, 
                    end = jends, 
                    chromo = GetChromo(jpeak),
                    genome = gen, name = "Data", type = "histogram")

plotTracks(list(itrack, gtrack, dtrack), collapseTranscripts="meta")
