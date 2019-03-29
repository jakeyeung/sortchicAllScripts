# Jake Yeung
# Date of Creation: 2019-03-28
# File: ~/projects/scchic/scripts/processing/lib/downstream_LDA_allele_specific.R
# Allele specific. Expects 40 chromosomes 

rm(list=ls())

timestart <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
# plotout <- paste0("/tmp/", jmark, "_LDA_bins_top_regions.pdf")

print(args)

jmark <- args[[1]]
inf <- args[[2]]
assertthat::assert_that(file.exists(inf))
plotout <- args[[3]]

nn <- as.numeric(args[[4]])
jmindist <- as.numeric(args[[5]])
sweep.params <- as.logical(args[[6]])

# sensible values
# nn=40
# jmindist=0.4
assertthat::assert_that(!is.na(nn))
assertthat::assert_that(!is.na(jmindist))
assertthat::assert_that(!is.na(jmindist))

pdf("/tmp/tmp_plots.pdf", useDingbats=FALSE)

# jmark <- "H3K4me1"
# jbin <- "TRUE"
# dirbase <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell"
# jdir <- paste0("lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin)
# if (jbin){
#   fname <- paste0("lda_out_meanfilt.BM-", jmark, ".CountThres0.K-5_10_15_20_25.Robj")
#   inf <- file.path(dirbase, jdir, fname)
# } else {
#   fname <- paste0("lda_out_meanfilt.BM-", jmark, ".CountThres0.K-5_15_25.Robj")  # I sweeped across fewer K's for jbin == FALSE, can add more later
#   inf <- file.path(dirbase, jdir, fname)
# }
# 

library(tidyr)  # gather

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


# Scripts for Allele-specific analysis ------------------------------------

SwapChromoFromCoord <- function(coord.old, chromo.hash, jsep = ":", outsep = "_"){
  # chr1:3000000-3100000 -> chr1_b6:2000000-31000000
  jsplit <- strsplit(coord.old, split = jsep)[[1]]
  chromo.old <- jsplit[[1]]
  coord.startend <- paste(jsplit[2:length(jsplit)], collapse = jsep)
  chromo.new <- chromo.hash[[chromo.old]]
  assertthat::assert_that(!is.null(chromo.new))
  coord.new <- paste(chromo.new, coord.startend, sep = jsep)
  return(coord.new)
}

StripSuffixFromCoord <- function(coord.full){
  # chrX_b6:48500000-48600000 -> chrX:48500000-48600000
  jsplit <- strsplit(coord.full, ":")[[1]]
  chromo <- jsplit[[1]]
  startend <- jsplit[[2]]
  chromo.new <- strsplit(chromo, "_")[[1]][[1]]
  return(paste0(chromo.new, ":", startend))
}

GetChromoSplitAllele <- function(x){
  xtmp <- GetChromo(x)
  return(strsplit(xtmp, "_")[[1]][[1]])
}



# Constants you can tweak -------------------------------------------------

# settings for H3K4me1
# the northern island with Sox6, Hbb, and Hba signal
# nearest neighbors settings for UMAP
nnvec <- seq(15, 51, 3)
mindistvec <- c(0.1, 0.2)
nnmindistvec <- apply(expand.grid(nnvec, mindistvec, stringsAsFactors=FALSE), 1, function(row) paste(row, collapse="_"))
names(nnmindistvec) <- nnmindistvec
print(nnmindistvec)


# other settings for UMAP same different marks
nnterms <- 15  # if you plot UMAP on the bins you use this variable
jmetric='euclidean' 
jseed=123


# LDA was run on binarized matrix or not. 
# I was thinking this binarized matrix would help reduce weird genomic regions with way too many reads. 
# Because we expect the count matrix to have only a few reads per bin per cell. 
# Can tweak this to TRUE or FALSE

top.thres <- 0.995


# Main code ---------------------------------------------------------------

print(paste("Loading LDA output", inf))

load(inf, v=T)

if (length(out.lda) > 1){
  out.lda <- ChooseBestLDA(out.lda)
  (kchoose <- out.lda@k)
} else {
  # choose first
  out.lda <- out.lda[[1]]
  kchoose <- out.lda@k
}

# Wrangle terms to Bcl6 vs Castaneous -------------------------------------

nchromos <- 20
chrs.raw <- seq(2 * nchromos)
chrs.new <- chrs.raw %% nchromos  # chromosome 0 is X chromosome
chrs.allele <- ceiling(chrs.raw / nchromos)
chrs.allele.named <- sapply(chrs.allele, function(x) ifelse(x == 1, "b6", "cast"))

chrs.dat <- data.frame(chromo.raw = paste("chr", chrs.raw, sep = ""),
                       chromo.num = as.character(chrs.raw), 
                       chromo.mod = paste("chr", chrs.new, sep = ""),
                       allele = chrs.allele.named, 
                       stringsAsFactors = FALSE)
chrs.dat$chromo.mod <- sapply(chrs.dat$chromo.mod, function(x) ifelse(x == "chr0", "chrX", x))
chrs.dat$chromo.full <- paste(chrs.dat$chromo.mod, chrs.dat$allele, sep = "_")

chromo.mod.hash <- hash::hash(chrs.dat$chromo.raw, chrs.dat$chromo.mod)
chromo.full.hash <- hash::hash(chrs.dat$chromo.raw, chrs.dat$chromo.full)
chromo.num.hash <- hash::hash(chrs.dat$chromo.num, chrs.dat$chromo.full)

# rename terms
terms.new <- sapply(out.lda@terms, SwapChromoFromCoord, chromo.full.hash, jsep = ":", outsep = "_")
terms.names.new <- sapply(names(out.lda@terms), SwapChromoFromCoord, chromo.num.hash, jsep = "_", outsep = "-")
names(terms.new) <- terms.names.new

# swap values of out.lda@terms
out.lda@terms <- terms.new


# Do topics ---------------------------------------------------------------



tm.result <- posterior(out.lda)

print(tail(colnames(tm.result$terms)))


topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

# for sweeping parameters
if (sweep.params){
    custom.settings.lst <- lapply(nnmindistvec, function(nn_mindist){
      nntmp <- as.numeric(strsplit(nn_mindist, "_")[[1]][[1]])
      mindisttmp <- as.numeric(strsplit(nn_mindist, "_")[[1]][[2]])
      settingstmp <- GetUmapSettings(nn=nntmp, jmetric=jmetric, jmindist=mindisttmp, seed = jseed)
      return(settingstmp)
    })

}
# for custom settings for real plot 
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

if (sweep.params){
    dat.umap.lst <- lapply(custom.settings.lst, function(csettings){
      dat.umap.tmp <- umap(topics.mat, config = csettings)
      rownames(dat.umap.tmp$layout) <- rownames(topics.mat)
      dat.umap.tmp <- data.frame(umap1 = dat.umap.tmp$layout[, 1], umap2 = dat.umap.tmp$layout[, 2])
      return(dat.umap.tmp)
    })
    names(dat.umap.lst) <- nnmindistvec
}

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)



# Plot dat umap -----------------------------------------------------------
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))


# Plot topic analysis -----------------------------------------------------


# Sort topics by interestingness ------------------------------------------

# analyze topic matrix across cells
topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>% 
  gather(key = "topic", value = "weight", -cell) %>%
  rowwise() %>%
  mutate(topic = as.numeric(substr(topic, 2, nchar(topic)))) %>%
  group_by(topic) %>%
  mutate(zscore = scale(weight, center = TRUE, scale = TRUE))

topics.sum <- topics.long %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(zscore < quantile(zscore, 0.99)) %>%
  mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)
print(topics.sum)

m1.cellweights <- ggplot(topics.long, aes(x = topic, y = zscore, group = topic)) + geom_boxplot() + theme_bw() + scale_x_continuous(breaks = seq(out.lda@k)) + ylab("Zscore Cell Weights") + ggtitle("For each cell, a datapoint for a bin weight in zscore. Many outliers suggest a group of high weight cells")

# plot distributions
m1.cellweights.hist <- ggplot(topics.long, aes(x = zscore, group = topic)) + geom_histogram(bins = 50) + facet_wrap(~topic) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("For each cell a dot for a bin weight (in zscore). Bimodal suggests a group of high weight cells")

barplot(height = topics.sum$entropy, names.arg = topics.sum$topic, xlab = "Topic", ylab = "Entropy Measure H", main = "Low H suggests topic with islands of high cell weights. High H evenly distributed weights", las = 2)

topics.ranked <- topics.sum$topic  # rank the topics by increasing entropy 

# Plot terms umap ---------------------------------------------------------


topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
top.regions <- unique(unlist(topic.regions))
terms.mat <- t(tm.result$terms)[top.regions, ]


# Peak analysis: annotate peaks propelry with the two alleles sepr --------


top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))
top.peaks$coord <- sapply(top.peaks$term, StripSuffixFromCoord)
regions <- data.frame(seqnames = sapply(colnames(tm.result$term), GetChromoSplitAllele),
                      start = sapply(colnames(tm.result$term), GetStart),
                      end = sapply(colnames(tm.result$term), GetEnd), 
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$term)

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)
top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))


# how top hits for specific topics
# Progenitor cellsare in Topic 12?
topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)
dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)

# plot an interesting region 

jtopic <- topics.ranked[[1]]  # toipc with lowest entropy 
jcol.rgb <- jcol.rgbs[[jtopic]]
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# make it pretty

m.interesting <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = paste0("`", jtopic, "`"))) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_continuous(paste0("Topic ", jtopic, "\nWeight")) + ggtitle("A topic with islands of cells with high weights")
print(m.interesting)


# translate beta to log fold change?
mat.norm <- t(tm.result$topics %*% tm.result$terms)  # this should give normalized signal, without the poisson noise?


cell.indx <- hash(rownames(dat.umap$knn$indexes), dat.umap$knn$indexes[, 1])
cell.indx.rev <- hash(dat.umap$knn$indexes[, 1], rownames(dat.umap$knn$indexes))
nr <- nrow(dat.umap$knn$indexes)
nc <- ncol(dat.umap$knn$indexes)
edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
colnames(edgelist) <- c("from", "to")
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


clstrs.orig <- as.character(sort(unique(as.numeric(dat.umap.long$louvain))))
# swap jclst with first element
clstrs.new <- clstrs.orig
remap.clstr <- hash(clstrs.orig, clstrs.new)
dat.umap.long$louvain <- sapply(as.character(dat.umap.long$louvain), function(x) remap.clstr[[x]])
dat.umap.long$louvain <- factor(as.character(dat.umap.long$louvain), levels = clstrs.orig)  # 1 to N
m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Spectral")
print(m.louvain)

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


gen <- "mm10"
chr <- "chr7"

mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")

dev.off()

pdf(plotout, useDingbats=FALSE)

if (sweep.params){
    # plot many UMAPs
    mlst <- lapply(nnmindistvec, function(x){
      dat.umap.tmp <- dat.umap.lst[[x]]
      jmain = x
      m <- ggplot(dat.umap.tmp, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmain)
      return(m)
    })
    lapply(mlst, print)
}

par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))

# order best topics
print(m1.cellweights)

# plot distributions
print(m1.cellweights.hist)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
barplot(height = topics.sum$entropy, names.arg = topics.sum$topic, xlab = "Topic", ylab = "Entropy Measure H", main = "Low H suggests topic with islands of high cell weights. High H evenly distributed weights", las = 2)


# show interesting topics

print(m.interesting)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
barplot(height = topics.sum$entropy, names.arg = topics.sum$topic, xlab = "Topic", ylab = "Entropy Measure H", main = "Topics are now ordered by increasing entropy.", las = 2)

# show top hits of betas
# for (i in seq(kchoose)){
for (i in topics.ranked){
  m.top <- subset(top.peaks.annotated, topic == i) %>% top_n(n=25, wt = beta) %>% 
    mutate(term = forcats::fct_reorder(term, dplyr::desc(beta))) %>%
    ggplot(aes(x = term, y = log10(beta), label = SYMBOL)) + 
    geom_point() + theme_bw(14) + geom_text_repel() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab("") + ylab("Log10 Bin Weight") + 
    ggtitle(paste("Top peak weights for topic:", i))
  print(m.top)
}

# print(m)

# plot topic weight
print(m.top)

# plot louvain clusters 
print(m.louvain)


dev.off()

print(Sys.time() - timestart)
