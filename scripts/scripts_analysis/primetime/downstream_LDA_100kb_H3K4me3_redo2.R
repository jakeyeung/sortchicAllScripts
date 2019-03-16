# Jake Yeung
# Date of Creation: 2019-03-05
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/downstream_LDA_100kb_H3K4me3_redo2.R
# H3K4me3


rm(list=ls())


jtime <- Sys.time() 

library(data.table)
library(tidyr)
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
library(forcats)
library(ggrepel)
library(biomaRt)

library(igraph)  # louvain

library(Gviz)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

source("scripts/Rfunctions/PlotFunctions.R")


# Load bulk RNAseq data ---------------------------------------------------

dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Constants you can tweak -------------------------------------------------

top.thres <- 0.995
inmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"
meanfilt <- 10
jbin <- "TRUE"; kstr <- "15_20_25_30_35"
# jbin <- "FALSE"; kstr <- "15_20_25_30"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
# jmark <- "H3K4me1"
jmark <- "H3K4me3"
# jmark <- "H3K9me3"

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

outmain <- "~/Dropbox/scCHiC_figs/FIG4_BM/primetime_plots"
plotout <- file.path(outmain, paste0(jmark, "_LDA_bins_top_regions.pdf"))

jtopic <- 1


topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

nn=35
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
jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)

print(sessionInfo())

# Plot dat umap -----------------------------------------------------------
# jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("pink", "darkred"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)


plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, asp = 0.75, cex = 0.2)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1), pty = "s")
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75, cex = 0.2)
}, jcol.rgbs, seq(kchoose))

# Plot terms umap ---------------------------------------------------------
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
top.regions <- unique(unlist(topic.regions))
terms.mat <- t(tm.result$terms)[top.regions, ]

# Uncomment below to plot the UMAP on the terms, can take a few more minutes

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
PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)

# how top hits for specific topics
# Progenitor cellsare in Topic 12?
topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)

# plot northern island
jcol.rgb <- jcol.rgbs[[jtopic]]
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# make it pretty

m.all <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2")) + geom_point(size = 0.8, alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
print(m.all)

# plot top hits for topic 12
print(subset(top.peaks.annotated, topic == 12), n = 50)
print(subset(top.peaks.annotated, topic == 10), n = 50)
print(subset(top.peaks.annotated, topic == 7), n = 50)

jgene <- c("Tgm2", "Mmrn1", "Trim48", "Pdzk1ip1", "Mllt3", "Mecom")
jgene <- c("Vldlr", "Uhrf1", "Rrm2", "Lig1", "Tipin")
jgene <- c("F2r", "Itga2b", "Zfp385a", "Zfpm1", "Plek", "Cd9", "Zeb2")

# translate beta to log fold change?
mat.norm <- t(tm.result$topics %*% tm.result$terms)  # this should give normalized signal, without the poisson noise?
# mat.norm <- mat.norm[top.regions, ]

# label using louvain clustering?

# do new knn

# dat.umap.fo

nn.louv=27
jmetric.louv='euclidean' 
jmindist.louv=0.4
jseed.louv=123

custom.settings.louv <- GetUmapSettings(nn=nn.louv, 
                                           jmetric=jmetric.louv, 
                                           jmindist=jmindist.louv, 
                                           seed = jseed.louv)

dat.umap.louv <- umap(topics.mat, config = custom.settings.louv)

dat.umap.louv.long <- data.frame(umap1 = dat.umap.louv$layout[, 1], umap2 = dat.umap.louv$layout[, 2], cell = rownames(dat.umap.louv$layout), 
                                 stringsAsFactors = FALSE)

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

dat.umap.long$louvain <- sapply(dat.umap.long$cell, function(x) clstr[[x]])

jclst <- 6
clstrs.orig <- as.character(sort(unique(as.numeric(dat.umap.long$louvain))))
# swap jclst with first element
clstrs.new <- clstrs.orig
clstrs.new[c(1, which(clstrs.new == jclst))] <- clstrs.new[c(which(clstrs.new == jclst), 1)]
remap.clstr <- hash(clstrs.orig, clstrs.new)


dat.umap.long$louvain <- sapply(as.character(dat.umap.long$louvain), function(x) remap.clstr[[x]])
dat.umap.long$louvain <- factor(as.character(dat.umap.long$louvain), levels = clstrs.orig)  # 1 to N
m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Spectral")
print(m.louvain)

# plot graph with edges?
# https://stackoverflow.com/questions/5364264/how-to-control-the-igraph-plot-layout-with-fixed-positions

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
coords <- layout.auto(g)
plot.igraph(igraph::simplify(g),
            layout = dat.umap$layout[V(g)$name, ],
            vertex.label = NA,
            edge.curved=FALSE,
            label = NA,
            edge.width = 0.5,
            vertex.size = 1)



# Merge cells and plot hits -----------------------------------------------


gen <- "mm10"
chr <- "chr7"

mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")

# get data for chromo 7
# redo for top 1000 peaks probably??

# jpeaks <- subset(top.peaks.annotated, rnk <= 100)$term
# jpeaks2 <- grep("chr7", rownames(mat.norm), value = TRUE)
# jpeaks <- c(jpeaks, jpeaks2)
jpeaks <- unique(top.peaks.annotated$term)
x <- as.data.frame(mat.norm[jpeaks, ])
x.long <- data.frame(exprs = unlist(x), cell = rep(colnames(x), each = nrow(x)), 
                     coord = rep(rownames(x), ncol(x)), stringsAsFactors = FALSE)
x.long$louvain.orig <- sapply(x.long$cell, function(x) clstr[[x]])
x.long$louvain <- x.long$louvain.orig
# x.long$louvain <- sapply(as.character(x.long$louvain.orig), function(x) remap.clstr[[x]])
x.long$exprs <- x.long$exprs * 10^6
x.long$coord <- sub("\\.", ":", x.long$coord)
x.long$coord <- sub("\\.", "-", x.long$coord)

# jpeaks <- grep("chr7", rownames(mat.norm), value = TRUE)
# x <- as.data.frame(mat.norm[jpeaks, ])
# x.long <- data.frame(exprs = unlist(x), cell = rep(colnames(x), each = nrow(x)), 
#                      coord = rep(rownames(x), ncol(x)), stringsAsFactors = FALSE)
# # x.long$louvain <- sapply(x.long$cell, function(x) clstr[[x]])
# x.long$louvain.orig <- sapply(x.long$cell, function(x) clstr[[x]])
# x.long$louvain <- sapply(as.character(x.long$louvain.orig), function(x) remap.clstr[[x]])
# x.long$exprs <- x.long$exprs * 10^6
# 
# # get data for chromo 11
# jpeaks11 <- grep("chr11", rownames(mat.norm), value = TRUE)
# x11 <- as.data.frame(mat.norm[jpeaks11, ])
# x.long11 <- data.frame(exprs = unlist(x11), cell = rep(colnames(x11), each = nrow(x11)), 
#                        coord = rep(rownames(x11), ncol(x11)), stringsAsFactors = FALSE)
# x.long11$louvain <- sapply(x.long11$cell, function(x) clstr[[x]])
# x.long11$exprs <- x.long11$exprs * 10^6



pdf(plotout, useDingbats = FALSE)

print(m.all)
print(m.louvain)

par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1), pty = "s")
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75, cex = 0.2)
}, jcol.rgbs, seq(kchoose))

par(mfrow=c(1, 1), mar=c(1,0.5,0.5,1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, asp = 0.75, cex = 0.2)
# plot topics
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75, cex = 0.2)
}, jcol.rgbs, seq(kchoose))

# show top hits of betas
topn <- 250
# show top hits of betas

# order topics before iterating

jtopics <- as.character(sort(unique(top.peaks.annotated$topic)))

# get average zscore 
top.genes.lst <- lapply(jtopics, function(jtopic){
  jtmp <- top.peaks.annotated %>% group_by(SYMBOL, topic) %>% filter(beta == max(beta)) %>% top_n(topn)
})

zscores.avg <- lapply(as.numeric(jtopics), function(jtopic){
  top.genes <- top.genes.lst[[jtopic]]$SYMBOL
  jsub <- subset(dat.long, Gene_Name %in% top.genes) %>%
    group_by(CellType) %>%
    summarise(zscore = mean(zscore)) %>%
    ungroup() %>%
    mutate(weight = exp(zscore) / sum(exp(zscore)))
  jsub$topic <- jtopic
  return(jsub)
}) %>%
  bind_rows()

# sort by range
zscores.ranked <- zscores.avg %>%
  group_by(topic) %>%
  summarise(diffrange = diff(range(zscore)),
            zscoremax = max(zscore)) %>%
  arrange(desc(diffrange))

zscores.ranked.bymax <- zscores.ranked %>%
  arrange(desc(zscoremax))

# jtopics.reordered <- zscores.ranked$topic
jtopics.reordered <- zscores.ranked.bymax$topic


for (i in as.character(jtopics.reordered)){
  m.top <- subset(top.peaks.annotated, topic == i) %>% top_n(n=25, wt = beta) %>% 
    mutate(term = forcats::fct_reorder(term, dplyr::desc(beta))) %>%
    ggplot(aes(x = term, y = log10(beta), label = SYMBOL)) + 
    geom_point() + theme_bw(14) + geom_text_repel() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab("") + ylab("Log10 Bin Weight") + 
    ggtitle(paste("Top peak weights for topic:", i))
  print(m.top)
  # plot RNAseq bulk gene expression
  top.genes <- unique(subset(top.peaks.annotated, topic == i)$SYMBOL[1:topn])
  jsub <- subset(dat.long, Gene_Name %in% top.genes)
  jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
  jlevels <- as.character(jsub.sorted.summarised$CellType)
  jsub$CellType <- factor(jsub$CellType, levels = jlevels)
  
  m.exprs <- ggplot(jsub, 
              aes(x = CellType , y = zscore)) + 
    geom_boxplot(outlier.shape = NA) +
    # geom_violin() +
    geom_jitter(width = 0.1, size = 0.5) +
    # geom_line() + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(paste("topic", i, "Top:", topn, "N Unique Genes", length(top.genes)))
  print(m.exprs)
  
}

# print(m)
# 
# # plot topic weight
# print(m.top)

# plot louvain clusters 
# print(m.louvain)

# plot peaks 


# # Sox6 region
# jgene <- "Sox6"
# jpeak <- subset(top.peaks.annotated, topic == jtopic & SYMBOL == jgene)$term[[1]]
# print(PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = FALSE, 
#                        return.plot.only = TRUE, usettings=custom.settings,
#                        gname = jgene))
# # for sox6 shift left a bit
# jstart <- subset(regions.annotated, region_coord == jpeak)$start - 6 * 10^5
# jend <- subset(regions.annotated, region_coord == jpeak)$end + 5 * 10^5
# PlotGTrack(x.long, jstart, jend, mart.obj, gen = "mm10", chr = "chr7", jheight = "auto")
# # PlotGTrack(x.long %>% mutate(louvain = ifelse(louvain == 1, "Eryth", "Others")), 
# #            jstart, jend, mart.obj, gen = "mm10", chr = "chr7", jheight = 1.5)
# 
# # Hbb region
# jgene <- "Hbb"
# jpeak <- subset(top.peaks.annotated, topic == jtopic & grepl(jgene, SYMBOL))$term[[1]]
# print(PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = FALSE, 
#                        return.plot.only = TRUE, usettings=custom.settings,
#                        gname = jgene))
# jstart <- subset(regions.annotated, region_coord == jpeak)$start - 5 * 10^5
# jend <- subset(regions.annotated, region_coord == jpeak)$end + 5 * 10^5
# 
# PlotGTrack(x.long, jstart, jend, mart.obj, gen = "mm10", chr = "chr7", jheight = "auto")
# # PlotGTrack(x.long %>% mutate(louvain = ifelse(louvain == 1, "Eryth", "Others")), 
# #            jstart, jend, mart.obj, gen = "mm10", chr = "chr7", jheight = 1.5)
# 
# # Hba region
# jgene <- "Hba"
# jchr <- "chr11"
# jpeak <- subset(top.peaks.annotated, topic == jtopic & grepl(jgene, SYMBOL))$term[[1]]
# print(PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = FALSE, 
#                        return.plot.only = TRUE, usettings=custom.settings,
#                        gname = jgene))
# jstart <- subset(regions.annotated, region_coord == jpeak)$start - 10 * 10^5
# jend <- subset(regions.annotated, region_coord == jpeak)$end + 10 * 10^5
# 
# PlotGTrack(x.long, jstart, jend, mart.obj, gen = "mm10", chr = jchr, jheight = "auto")
# 
# # PlotGTrack(x.long %>% mutate(louvain = ifelse(louvain == 1, "Eryth", "Others")),
# #            jstart, jend, mart.obj, gen = "mm10", chr = jchr, jheight = 1.5)
# 
# # Do top peak for every region
# 

for (jtopic in jtopics.reordered){
  jpeak <- subset(top.peaks.annotated, topic == jtopic)$term[[1]]
  print(PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = FALSE, 
                         return.plot.only = TRUE, usettings=custom.settings,
                         gname = jgene))
  jstart <- subset(regions.annotated, region_coord == jpeak)$start - 5 * 10^5
  jend <- subset(regions.annotated, region_coord == jpeak)$end + 5 * 10^5
  jchromo <- subset(regions.annotated, region_coord == jpeak)$seqnames
  
  if (jchromo == "chr20" | jchromo == "chr21"){
    print("Skipping sex chromos")
    next 
 }
  
  PlotGTrack(x.long, jstart, jend, mart.obj, gen = "mm10", chr = jchromo, jheight = "auto")
  # PlotGTrack(x.long %>% mutate(louvain = ifelse(louvain == 1, "Eryth", "Others")), 
  #            jstart, jend, mart.obj, gen = "mm10", chr = jchromo, jheight = 1.5)
}
dev.off()
print(Sys.time() - jtime)

