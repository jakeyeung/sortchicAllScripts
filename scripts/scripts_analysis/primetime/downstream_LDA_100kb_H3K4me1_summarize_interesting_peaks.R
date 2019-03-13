# Jake Yeung
# Date of Creation: 2019-03-05
# Jake Yeung
# Date of Creation: 2019-03-05
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/downstream_LDA_100kb_H3K4me1_summarize_interesting_peaks.R
# Show interesting peaks related to specific cell types 

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

GetTopGenes <- function(dat.long, jcelltype, jtopic, topn = 20){
  jsub <- dat.long %>% 
    rowwise() %>%
    mutate(CellType2 = ifelse(CellType == "nucleate_erythrocyte", TRUE, FALSE)) %>% 
    group_by(CellType2, Gene_Name) %>%
    summarise(zscore = mean(logFPKM)) %>%
    group_by(Gene_Name) %>% 
    summarise(zscore.diff = diff(range(zscore)))
  # jsub.hash <- hash(jsub$Gene_Name, jsub$zscore.diff)
  # get top genes
  jsub.top <- top.peaks.annotated %>% 
    filter(topic == jtopic & rnk <= 100) %>%
    group_by(topic, SYMBOL) %>%
    filter(beta == max(beta)) %>%
    arrange(rnk)
  jsub.top <- left_join(jsub.top, jsub, by = c("SYMBOL" = "Gene_Name")) %>% 
    arrange(desc(zscore.diff))
  jsub.sub <- jsub.top[1:topn, ]
  # genes <- jsub.top$SYMBOL[1:topn]
  return(jsub.sub)
}

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
jmark <- "H3K4me1"
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
plotout <- file.path(outmain, paste0(jmark, "_LDA_summarized_bins.pdf"))

jtopic <- 1


topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

nn=40
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


# Cherry pick regions  ----------------------------------------------------

head(top.peaks.annotated)

jtopn <- 5
# erythrypoiesis peaks, topic 1
jtopic <- 1
jcelltype <- "nucleate_erythryocyte"
genes.eryth <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)

# granulocyte peaks, topic 21 or topic 13?
# jtopic <- 21
jtopic <- 13
jcelltype <- "granulocyte"
genes.gran <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
genes.gran2 <- GetTopGenes(dat.long, jcelltype, 13, topn = jtopn)

# B-cells
# jtopic <- 22
jtopic <- 11
jcelltype <- "lymphocyte_of_B_lineage"
genes.b <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
# genes.bcells <- c("Irf4", "Irf8", "Il7r", "Ppp1r16b", "Bach2")

# T-cells
jtopic <- 14
jcelltype <- "T_cell"
# 14 is T-cells
genes.t <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `14`)) + geom_point() 
# interesting genes: Ly6c2 by https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4039991/ Teichmann

# Natural killer cells
jtopic <- 27
jcelltype <- "natural_killer_cell"
genes.nk <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `27`)) + geom_point() 

# Megakaryocytes
jtopic <- 15
jcelltype <- "megakaryocyte"
genes.mega <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `15`)) + geom_point() 

# Monocytes
# jtopic <- 13
jtopic <- 5
jcelltype <- "monocyte"
genes.mono <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `5`)) + geom_point() 

# Dendritic cells
jtopic <- 9
jcelltype <- "dendritic_cell"
genes.den <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `9`)) + geom_point() 

# Stem cells??
jtopic <- 12
jcelltype <- "Kit_and_Sca1−positive_hematopoietic_stem_cell"
genes.stem <- GetTopGenes(dat.long, jcelltype, jtopic, topn = jtopn)
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `12`)) + geom_point() 


# Heatmap of peaks?  ------------------------------------------------------

# heatmap of peak weights for each topic 
peaks <- c(genes.eryth$term, genes.gran$term, genes.b$term, genes.t$term, genes.nk$term, genes.mega$term, genes.mono$term, genes.den$term, genes.stem$term)
# add Hbb
peaks.more <- subset(top.peaks.annotated, grepl("Hbb|Hba", SYMBOL))[1:2, ]$term
peaks <- c(peaks, peaks.more)

peaks.i <- which(colnames(tm.result$terms) %in% peaks)
topics <- c(1, 13, 11, 14, 27, 15, 5, 9, 12)
topics.i <- which(rownames(tm.result$terms) %in% topics)
terms.sub <- tm.result$terms[topics.i, peaks.i]

# rename peaks
peaks.hash <- hash(top.peaks.annotated$term, paste(top.peaks.annotated$term, top.peaks.annotated$SYMBOL, sep = ";"))

colnames(terms.sub) <- sapply(colnames(terms.sub), function(x) peaks.hash[[x]])

# plot heatmap
library(heatmap3)
library(gplots)
library(made4)

pdf("~/Dropbox/scCHiC_figs/FIG4_BM/primetime_plots/topics_and_heatmap.pdf", useDingbats = FALSE)
# plot topics

for (topic in topics){
  topic <- as.character(topic)
  dat.tmp <- dat.umap.long %>% dplyr::select(c(umap1, umap2, topic))
  dat.tmp[[topic]] <- CapQuantile(dat.tmp[[topic]], cap.quantile = 0.99)
  midpt <- min(dat.tmp[[topic]]) + (max(dat.tmp[[topic]]) - min(dat.tmp[[topic]])) / 2
  # compress 
  # print(midpt)
  dat.tmp <- RankOrder(dat.tmp, topic)
  m.top <- ggplot(dat.tmp, aes_string(x = "umap1", y = "umap2", color = paste0("`", topic, "`")), order = orderrank) + geom_point() +
    scale_color_gradient2(low = "gray95", mid = "gray50", high = "darkred", midpoint = midpt)
  print(m.top)
}

# pl t hits

heatmap3(t(terms.sub), scale = "row", margins = c(5, 8), cexRow = 0.25)
dev.off()
