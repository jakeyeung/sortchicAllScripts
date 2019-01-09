# Jake Yeung
# Date of Creation: 2019-01-08
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_analysis_K562_cell_cycle.R
# Look at cell cycle

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(destiny)
library(hash)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")



# Functions ---------------------------------------------------------------

ExtractCellCycle <- function(x){
  
}


# Load --------------------------------------------------------------------

# inf <- "/private/tmp/K562/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE/lda_out_meanfilt.PZ-K562-H3K4me3.CountThres0.K-10_15_20_25.Robj"
# inf <- "/private/tmp/K562/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE/lda_out_meanfilt.PZ-K562-H3K4me3.CountThres0.K-10_15_20_25.Robj"
# inf <- "/private/tmp/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-5_7_10_12_15_20_25_30.Robj"
# inf <- "/private/tmp/lda_out_meanfilt.PZ-BM-H3K4me3.CountThres0.K-5_7_10_12_15_20_25_30.Robj"
# inf <- "/private/tmp/lda_out_meanfilt.PZ-BM-H3K9me3.CountThres0.K-5_7_10_12_15_20_25_30.Robj"
# inf <- "/private/tmp/lda_out_meanfilt.PZ-BM-H3K27me3.CountThres0.K-5_7_10_12_15_20_25_30.Robj"
# inf <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K4me3.CountThres0.K-10_15_20_25.Robj"
inf <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K27me3.CountThres0.K-10_15_20_25.Robj"

load(inf, v=T)



# get count.mat
# count.mat <- GetCountMatFromLDA(out.lda[[1]])

print(1 - Matrix::nnzero(count.mat) / length(count.mat))

out.lda.lst <- out.lda

# we did multicore, so which one to choose?
Kvec <- sapply(out.lda.lst, function(x) x@k)

# plot best fit?
# result <- CalculateMetrics(out.lda.lst, Kvec)

# print(result)

# FindTopicsNumber_plot(result)

best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]

# plot loglikelihood
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood), 'o')

# plot likelihood 

# kchoose <- best.K
kchoose <- 20
out.lda <- out.lda.lst[[which(Kvec == kchoose)]]


# Get cell cycle info -----------------------------------------------------

cellnames <- unname(out.lda@documents)
cellcycles <- sapply(cellnames, function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE)
print(unique(cellcycles))
colhash <- hash(c("G1", "G2M", "S"), c("red", "blue", "orange"))
jcol.phase <- sapply(cellcycles, function(x) colhash[[x]], USE.NAMES = FALSE)

tmResult <- posterior(out.lda)

nn <- 5
# jmetric <- 'pearson'
# jmetric <- 'cosine'
jmetric <- 'euclidean'
jmindist <- 0.0001
custom.settings <- umap.defaults
custom.settings$n_neighbors <- nn
custom.settings$metric <- jmetric
custom.settings$min_dist <- jmindist

dat.umap <- umap(tmResult$topics, config = custom.settings)
rownames(dat.umap$layout) <- rownames(tmResult$topics)

jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# color by Erdr1 reads
jcounts <- Matrix::colSums(count.mat)
jcol.counts <- ColorsByCounts(jcounts, nbreaks = 100)

par(mfrow=c(1,1), mar=c(1,1,1,1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col = jcol.phase)

dat.pca <- prcomp(tmResult$topics, center = TRUE, scale. = TRUE)

plot(dat.pca$x[, 1], dat.pca$x[, 2], col = jcol.phase, pch = 20)
# color by gamma

jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma)

nb.col <- 5
nb.row <- ceiling(best.K / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))


# Destiny ------------------------------------------------------------------

dm.out <- destiny::DiffusionMap(data = tmResult$topics)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dm.out, 1:2, col = jcol.phase)  # colors from beta matrix
# plot(dm.out, 1:3, col = jcol.phase)  # colors from beta matrix


(x <- terms(out.lda, 20))

# Hox genes?
hoxa <- "chr16:521"
hoxb <- "chr11:962"
hoxc <- "chr15:1030"
hoxc2 <- "chr15:1029"
hoxd <- "chr2:747"
hbb <- "chr7:1038"

greplst <- list(hoxa, hoxb, hoxc, hoxc2, hoxd, hbb)
x.match <- lapply(greplst, function(grepstr) grep(grepstr, x, value = TRUE))
print(x.match)


# Do GREAT analysis -------------------------------------------------------

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

library(ChIPseeker)
library(rGREAT)
library(JFuncs)
library(gwascat)

# library(rtracklayer)

top.thres <- 0.98
# need to assign cutoff for each peak for each topic
topic.regions <- lapply(seq(best.K), function(clst){
  return(SelectTopRegions(tmResult$terms[clst, ], colnames(tmResult$terms), method = "thres", method.val = top.thres))
})
regions <- data.frame(seqnames = sapply(colnames(tmResult$terms), GetChromo),
                      start = sapply(colnames(tmResult$terms), GetStart),
                      end = sapply(colnames(tmResult$terms), GetEnd))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))

# annotate after HG19
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                                annoDb='org.Hs.eg.db'))
rownames(regions.annotated) <- regions.annotated$region_coord

ncores <- 1

# liftover regions from Hg38 to Hg19
# https://bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- rtracklayer::import.chain(path)

seqlevelsStyle(gr.in) = "UCSC"
regions.range.19 = unlist(rtracklayer::liftOver(regions.range, ch))
regions.range.19$hg38peak <- names(regions.range.19)

# annotate regions with hg19
# annotate after HG19
i <- 1
# gr.in <- regions.range.19[topic.regions[[i]], ]
# gr.in <- subset(regions.range.19, hg38peak %in% topic.regions[[i]])
regions.annotated.19 <- as.data.frame(annotatePeak(unname(regions.range.19), 
                                                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                                annoDb='org.Hs.eg.db'))
# regions.annotated.19$hg38peak <- names(regions.range.19)
# 
# print(i)
# print(head(gr.in))
gr.in <- subset(regions.range.19, hg38peak %in% topic.regions[[i]])
out.great <- submitGreatJob(unname(gr.in), species="hg19", request_interval = 10)
out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great),
                              request_interval = 300)


if (ncores == 1){
  out.great.lst <- lapply(seq(best.K), function(i){
    gr.in <- regions.range[topic.regions[[i]], ]
    print(i)
    print(head(gr.in))
    out.great <- submitGreatJob(gr.in, species="hg19", request_interval = 700)
    return(out.great)
  })
  out.tb.lst <- lapply(out.great.lst, function(out.great){
    out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great),
                                  request_interval = 300)
    return(out.tb)
  })
} else {
  print(paste("Running great multicore", ncores))
  out.great.lst <- mclapply(seq(best.K), function(i){
    gr.in <- regions.range[topic.regions[[i]], ]
    out.great <- submitGreatJob(gr.in, species="hg19", request_interval = 700)
    return(out.great)
  }, mc.cores = ncores)
  out.tb.lst <- mclapply(out.great.lst, function(out.great){
    out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great), 
                                  request_interval = 300)
    return(out.tb)
  }, mc.cores = ncores)
}
