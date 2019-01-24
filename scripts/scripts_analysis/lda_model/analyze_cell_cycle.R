# Jake Yeung
# Date of Creation: 2019-01-15
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/analyze_cell_cycle.R
# Analyze cell cycle from count matrix

rm(list=ls())

library(dplyr)
library(ggplot2)
library(topicmodels)
library(cisTopic)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(hash)
library(umap)

source("scripts/Rfunctions/ParseStrings.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")



# Load data ---------------------------------------------------------------

# jchip <- "H3K4me3"
jchip <- "H3K4me1"
# jchip <- "H3K27me3"
dirmain <- "/tmp/count_mat_K562_round2"
inf <- file.path(dirmain, paste0("count_mats.hiddenDomains.1000/PZ-K562-", jchip, ".merged.hiddenDomains.NoCountThres.Robj"))

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

count.mat <- count.dat$counts

print(dim(count.mat))

dat.meanvar <- data.frame(Sum = Matrix::rowSums(count.mat), 
                          Mean = Matrix::rowMeans(count.mat),
                          Var = apply(count.mat, 1, var),
                          peak = rownames(count.mat),
                          stringsAsFactors=FALSE)
dat.meanvar <- dat.meanvar %>%
  rowwise() %>%
  mutate(CV = sqrt(Var) / Mean,
         peaksize = GetPeakSize(peak))

p1 <- ggplot(dat.meanvar, aes(x = log10(peaksize))) + geom_histogram(bins=25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


p2 <- ggplot(dat.meanvar, aes(x = log10(Mean), y = log10(CV), size = peaksize)) + geom_point(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = -0.5)

p3 <- ggplot(dat.meanvar, aes(x = log10(peaksize), y = Sum)) + geom_point(alpha = 0.1) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

print(p1)
print(p2)
print(p3)

# Limit to cell cycle genes -----------------------------------------------

glst <- "/Users/yeung/projects/scChiC/data/genelists/cyclebaseGeneList"
cc.genes.dat <- data.table::fread(glst)
cc.genes <- cc.genes.dat$Genename

# get regions
terms <- rownames(count.mat)
regions <- data.frame(seqnames = sapply(terms, GetChromo),
                      start = sapply(terms, GetStart),
                      end = sapply(terms, GetEnd))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))

# chainpath <- "~/scChiC/data/chainfiles/hg38ToHg19.over.chain"
chainpath <- "data/chainfiles/hg38ToHg19.over.chain"
ch <- rtracklayer::import.chain(chainpath)
seqlevelsStyle(regions.range) = "UCSC"
regions.range.19 = unlist(rtracklayer::liftOver(regions.range, ch))
regions.range.19$hg38peak <- names(regions.range.19)

regions.annotated.19 <- as.data.frame(annotatePeak(unname(regions.range.19),
                                                   TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                                                   annoDb='org.Hs.eg.db'))

regions.filt <- subset(regions.annotated.19, SYMBOL %in% cc.genes & abs(distanceToTSS) < 100000)

peaks.keep <- regions.filt$hg38peak

peaks.keep.i <- which(rownames(count.mat) %in% peaks.keep)

# check this
bad.cell <- "S_AH2VV5BGX9_S9.filtered.sorted.TAGTACGT.sorted.bam-K562-H3K4me3-NA-AH2VV5BGX9-TAGTACGT"

# filter for CC genes
count.mat <- count.mat[peaks.keep.i, ]


# Binarize and go? --------------------------------------------------------

meanmax <- 1
cellmin <- 100
cellmax <- 50000

# suspicious peaks
bad.peaks <- dat.meanvar %>%
  filter(Mean > meanmax)
print(head(bad.peaks$peak))
print(paste("There are", nrow(bad.peaks), "peaks with more than", "counts. Removing them..."))

# filter them out before running count.mat
print("Dimensions before filtering peaks...")
print(dim(count.mat))

count.mat <- count.mat[which(!rownames(count.mat) %in% bad.peaks$peak), ]
print("Dimensions after filtering peaks...")
print(dim(count.mat))

# remove peaks in chrM
M.peaks <- grep("chrM|chrY|chrX", dat.meanvar$peak, value=TRUE)
# M.peaks <- grep("chrM", dat.meanvar$peak, value=TRUE)

count.mat <- count.mat[which(!rownames(count.mat) %in% M.peaks), ]
print("Dimensions after filtering peaks X, Y, M chromos")
print(dim(count.mat))


# Remove cells with zero entries
print("Dimensions before filtering cells")
print(dim(count.mat))
count.mat <- count.mat[, which(Matrix::colSums(count.mat) > cellmin)]
count.mat <- count.mat[, which(Matrix::colSums(count.mat) < cellmax)]
print("Dimensions after filtering cells...")
print(dim(count.mat))


# Binarize matrix and run -------------------------------------------------

# count.mat <- BinarizeMatrix(count.mat)

# fix column names
colnames(count.mat) <- unname(colnames(count.mat))

# add meta data
cellcycles <- unname(sapply(colnames(count.mat), function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE))
print(unique(cellcycles))

metadat <- data.frame(cname=colnames(count.mat), cycle=cellcycles, stringsAsFactors = FALSE)
rownames(metadat) <- metadat$cname

nclst <- 5
print("Running single LDA for topics:")
print(nclst)

# # maybe try to run cistopics?
# cisTopicObject <- createcisTopicObject(count.mat, project.name=paste0("K562", jchip, "cellcycle"))
# cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadat)
# 
# system.time(
#   cisTopicObject <- runModels(cisTopicObject, topic=c(10, 15, 20, 25), seed=987, nCores=4, burnin = 120, iterations = 150, addModels=FALSE)
# )
# 
# cisTopicObject <- selectModel(cisTopicObject)
# 
# cisTopicObject <- runUmap(cisTopicObject, target='cell', n_components=3)
# cisTopicObject <- runPCA(cisTopicObject, target='cell')
# cisTopicObject <- runtSNE(cisTopicObject, target='cell', dim=3)
# 
# par(mfrow=c(1,1))
# 
# jcols <- c("blue", "red", "brown")
# jcols <- list(G1='blue', G2M='red', S='brown')

# plotFeatures(cisTopicObject, method='PCA', target='cell', topic_contr=NULL,
#              colorBy=c('cycle'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE,
#              col.low='black', col.mid='orange', col.high='red', intervals=20)
# plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, 
#              colorBy=c('cycle'), cex.legend = 0.8, factor.max=.75, dim=3, legend=TRUE, intervals=20)
# plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, 
#              colorBy=c('cycle'), cex.legend = 0.8, factor.max=.75, dim=3, legend=TRUE, 
#              col.low='black', col.mid='orange', col.high='red', intervals=20)

system.time(
  out.lda <- LDA(x = t(count.mat), k = nclst, method = "Gibbs", control=list(seed=0))
)

kchoose <- out.lda@k
tm.result <- posterior(out.lda)

cellnames <- unname(out.lda@documents)
cellcycles <- sapply(cellnames, function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE)
print(unique(cellcycles))

phasenames <- c("G1", "G2M", "S")
phasecols <- c("red", "blue", "black")
colhash <- hash(phasenames, phasecols)
jcol.phase <- sapply(cellcycles, function(x) colhash[[x]], USE.NAMES = FALSE)

nn <- 4
nn.terms <- 15
# jmetric <- 'pearson'
# jmetric <- 'cosine'
jmetric <- 'euclidean'
jmindist <- 0.5
custom.settings <- GetUmapSettings(nn, jmetric, jmindist)
custom.settings.terms <- GetUmapSettings(nn.terms, jmetric, jmindist)

dat.umap <- umap(tm.result$topics, config = custom.settings)
rownames(dat.umap$layout) <- rownames(tm.result$topics)

jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col = jcol.phase)
legend(2, -4, legend=phasenames, col=phasecols, pch = 20)

dat.pca <- prcomp(tm.result$topics, center = TRUE, scale. = TRUE)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dat.pca$x[, 1], dat.pca$x[, 2], col = jcol.phase, pch = 20)


