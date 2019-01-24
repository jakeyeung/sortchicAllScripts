# Jake Yeung
# Date of Creation: 2019-01-21
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/analyze_cell_cycle_round2_downstream.R
# 

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
library(plotly)

source("scripts/Rfunctions/ParseStrings.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/GetMetaCellHash.R")



# Get cell cycle metadata -------------------------------------------------

CellCycleHash <- function(n.rows = 16, n.cols = 24){
  # first 1/3 of plate is G1, 2/d is S, 3/3 is G2
  nsplits <- 3
  width <- n.cols / nsplits
  cellids <- rep(NA, n.rows * n.cols)
  for (i in 1:n.rows){
    indx.G1 <- (1 : width) + (i - 1) * (width * nsplits)
    indx.S <- ((width+1) : (width*2)) + (i - 1) * (width * nsplits)
    indx.G2 <- ((width*2+1) : (width*3)) + (i - 1) * (width * nsplits) 
    cellids[indx.G1] <- "G1"
    cellids[indx.S] <- "S"
    cellids[indx.G2] <- "G2"
  }
  # return as hash
  jnames <- paste("cell", seq(length(cellids)), sep = "")
  cellhash <- hash(jnames, cellids)
  return(cellhash)
}

cchash <- CellCycleHash()



# Load data ---------------------------------------------------------------

# jchip <- "H3K4me3"
# jchip <- "H3K4me1"
jchip <- "H3K27me3"

dirmain <- "/tmp/count_mat_K562_round2_LDA_output"
inf <- file.path(dirmain, paste0("lda_outputs.hiddenDomains.meanfilt_1.cellmin_100.cellmax_500000.binarize.FALSE/lda_out_meanfilt.PZ-K562-", jchip, ".CountThres0.K-10_15_20_25.Robj"))
inf.GREAT <- file.path(dirmain, paste0("lda_outputs.hiddenDomains.meanfilt_1.cellmin_100.cellmax_500000.binarize.FALSE/downstream/lda_out_meanfilt.PZ-K562-", jchip, ".CountThres0.K-10_15_20_25.GREAT.Robj"))

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.GREAT))

load(inf, v=T)
load(inf.GREAT, v=T)



# handle round2 colnames
barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))

cnames.old <- unname(colnames(count.mat))
repl <- sapply(cnames.old, function(x) strsplit(x, "-")[[1]][[4]])
bcs <- sapply(cnames.old, function(x) strsplit(x, "-")[[1]][[6]])
cellid <- sapply(bcs, function(x) cellhash.bc[[x]])
cellcycle <- sapply(cellid, function(x) cchash[[x]])
name2phase <- hash(cnames.old, cellcycle)
# Get plate and barcode

unique(sapply(cnames.old, function(x) strsplit(x, "-")[[1]][[4]]))

colnames(count.mat) <- unname(colnames(count.mat))

# add meta data
cellcycles <- sapply(unname(colnames(count.mat)), function(x) name2phase[[x]])
print(unique(cellcycles))

metadat <- data.frame(cname=colnames(count.mat), cycle=cellcycles, stringsAsFactors = FALSE)
rownames(metadat) <- metadat$cname



kchoose <- out.lda@k
tm.result <- posterior(out.lda)

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


phasenames <- c("G1", "G2", "S")
phasecols <- c("red", "blue", "cyan")
colhash <- hash(phasenames, phasecols)
jcol.phase <- sapply(cellcycles, function(x) colhash[[x]], USE.NAMES = FALSE)

pdf(paste0("~/Dropbox/scCHiC_figs/FIG2_K562-G1/cell_cycle/", jchip, ".cellcycle.pdf"))

dat.pca <- prcomp(tm.result$topics, center = TRUE, scale. = TRUE)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dat.pca$x[, 1], dat.pca$x[, 2], pch = 20, col = jcol.phase, main = jchip)
legend(2, -4, legend=phasenames, col=phasecols, pch = 20)
plot(dat.pca$x[, 2], dat.pca$x[, 3], pch = 20, col = jcol.phase, main = jchip)
legend(2, -4, legend=phasenames, col=phasecols, pch = 20)

out.pca <- as.data.frame(dat.pca$x)
out.pca$phase <- sapply(rownames(out.pca), function(x){
  x <- strsplit(x, "-")[[1]][[6]]
  cchash[[cellhash.bc[[x]]]]
})
m <- plotly::plot_ly(out.pca, x=~PC1, y=~PC2, z=~PC3, color=~phase, colors = c('red', 'blue', 'cyan'), 
                     size=1) %>% layout(title=jchip)
print(m)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col=jcol.phase)
legend(2, -4, legend=phasenames, col=phasecols, pch = 20)

print(m)


mat.norm <- t(tm.result$topics %*% tm.result$terms)
# Plot terms --------------------------------------------------------------

topics.mat <- tm.result$topics

jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.pca$x[, 1], dat.pca$x[, 2], pch = 20, col = jcol.rgb, main = jtopic)
  # plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))


library(tidytext)
ap_topics <- tidy(out.lda, matrix = "beta") %>% arrange(desc(beta))
print(subset(ap_topics, topic == 10))
print(subset(ap_topics, topic == 15))

jsub <- subset(ap_topics, topic %in% c(10, 15) & grepl("chr1", term)) %>% 
  group_by(term) %>%
  filter(beta == max(beta)) %>%
  arrange(term)

# topic 10
jpeak <- "chr9:113293000-113336000"  # WDR31
jpeak <- "chr1:9617000-9655000"  # PIK3CD

# topic 15
jpeak <- "chr16:83498000-83546000"  # CDH13

jpeaks <- c("chr9:113293000-113336000", "chr1:9617000-9655000", "chr16:83498000-83546000")

jpeak <- "chr1:100528000-100554000"

for (jpeak in jpeaks){
  # plot topic 10 hits 
  peaks.keep <- jpeak
  row.i <- which(rownames(mat.norm) %in% peaks.keep)
  if (length(row.i) == 0){
    print(paste("Skipping peak", peaks.keep))
    next()
  }
  
  jcounts.norm <- mat.norm[row.i, ]
  print(range(jcounts.norm))
  
  jcol.counts <- ColorsByCounts(jcounts.norm, nbreaks = 100, colvec = c("lightblue", "darkblue"))
  
  # plot topics soft clustering weights
  topics.mat <- tm.result$topics
  jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
  nb.col <- 5
  nb.row <- ceiling(kchoose / nb.col)
  
  
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  plot(dat.pca$x[, 1], dat.pca$x[, 2], pch = 20, col = jcol.counts, main = paste("Norm counts", jpeak))
  
}


dev.off()





