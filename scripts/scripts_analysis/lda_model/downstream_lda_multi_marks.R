# Jake Yeung
# Date of Creation: 2019-01-06
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_lda_multi_marks.R
# Analyze downstream of multiple marks

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)


source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")


pcut <- 0.3
jdist <- 1000
jchip <- "H3K4me1"
jchip <- "H3K9me3"
jchip <- "H3K4me3"
jchip <- "H3K27me3"
load(paste0('/private/tmp/lda_output.systematic.', pcut, '.', jdist, '/lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-5_7_10_12_15_20_25_30.Robj'), v=T)

# load('/private/tmp/lda_output.systematic.0.3.1000/lda_out_meanfilt.PZ-BM-H3K27me3.CountThres0.K-5_7_10_12_15_20_25_30.Robj', v=T)
# load('/private/tmp/lda_output.systematic.0.3.1000/lda_out_meanfilt.PZ-BM-H3K9me3.CountThres0.K-5_7_10_12_15_20_25_30.Robj', v=T)
# load('/private/tmp/lda_output.systematic.0.3.1000/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-5_7_10_12_15_20_25_30.Robj', v=T)
# load('/private/tmp/lda_output.systematic.0.3.1000/lda_out_meanfilt.PZ-BM-H3K4me3.CountThres0.K-5_7_10_12_15_20_25_30.Robj', v=T)


# Plot outputs ------------------------------------------------------------

# rebuild count matrix

out.lda.lst <- out.lda

# we did multicore, so which one to choose?
Kvec <- sapply(out.lda.lst, function(x) x@k)

# plot best fit?
# result <- CalculateMetrics(out.lda.lst, Kvec)

# print(result)

# FindTopicsNumber_plot(result)

best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]

# plot loglikelihood
par(mfrow=c(1,1), mar=c(1,1,1,1))
plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood))

# plot likelihood 

kchoose <- best.K
out.lda <- out.lda.lst[[which(Kvec == kchoose)]]

tmResult <- posterior(out.lda)

nn <- 15
# jmetric <- 'pearson2'
jmetric <- 'euclidean'
# jmetric <- 'cosine'
jmindist <- 0.001
custom.settings <- umap.defaults
custom.settings$n_neighbors <- nn
custom.settings$metric <- jmetric
custom.settings$min_dist <- jmindist

dat.umap <- umap(tmResult$topics, config = custom.settings)
rownames(dat.umap$layout) <- rownames(tmResult$topics)

jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# color by Erdr1 reads
(jpeaks1 <- grep("chrY:90[7,8,9][0-9]{5}-", rownames(count.mat), value = TRUE))  # Erdr1 on Y
(jpeaks2 <- grep("chrX:1699{6}-", rownames(count.mat), value = TRUE))  # Erdr1 on X
(jpeaks3 <- grep("chrX:170[0-9]{6}-", rownames(count.mat), value = TRUE))  # Erdr1 on X
jpeaks.all <- c(jpeaks1, jpeaks2, jpeaks3)
# jpeaks.all <- c(jpeaks1)
jpeaks.i <- which(rownames(count.mat) %in% jpeaks.all)
jcounts <- Matrix::colSums(count.mat[jpeaks.i, ])
jcol.counts <- ColorsByCounts(jcounts, nbreaks = 100)

par(mfrow=c(1,1), mar=c(1,1,1,1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col = jcol.counts)


# color by gamma

jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma)

nb.col <- 5
nb.row <- ceiling(best.K / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))

(x <- terms(out.lda, 10))

# color by top peaks
jpeaks <- terms(out.lda, 20)[1:20, 20]
# (jpeaks <- grep("chrY:90", rownames(count.mat), value = TRUE))  # Erdr1 genes
jpeaks.i <- which(rownames(count.mat) %in% jpeaks)
jcounts <- Matrix::colSums(count.mat[jpeaks.i, ])
jcol.counts <- ColorsByCounts(jcounts, nbreaks = 100)
par(mfrow=c(1,1), mar=c(1,1,1,1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col = jcol.counts)


# Hox genes?
hoxa <- "chr16:52"
hoxb <- "chr11:96"
hoxc <- "chr15:102"
hoxc2 <- "chr15:103"
hoxd <- "chr2:74"

greplst <- list(hoxa, hoxb, hoxc, hoxc2, hoxd)

x.match <- lapply(greplst, function(grepstr) grep(grepstr, x, value = TRUE))

print(x.match)


# Problem peaks -----------------------------------------------------------

bad.peak <- "chr8:128084965-128089512"
peak.i <- which(rownames(count.mat) == bad.peak)

slice <- count.mat[peak.i, ]

# plot average across peaks
avgs <- sort(Matrix::rowMeans(count.mat), decreasing = TRUE)

par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(density(log10(avgs)))
plot(density(avgs))


# Run GREAT --------------------------------------------------------------

# need to assign cutoff for each peak for each topic
top.thres <- 0.98
topic.regions <- lapply(seq(best.K), function(clst){
  return(SelectTopRegions(tmResult$terms[clst, ], colnames(tmResult$terms), method = "thres", method.val = top.thres))
})

regions <- data.frame(seqnames = sapply(colnames(tmResult$terms), GetChromo),
                      start = sapply(colnames(tmResult$terms), GetStart),
                      end = sapply(colnames(tmResult$terms), GetEnd))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
rownames(regions.annotated) <- regions.annotated$region_coord

# takes ~10 minutes

out.great.lst <- lapply(1, function(i){
  gr.in <- regions.range[topic.regions[[i]], ]
  out.great <- submitGreatJob(gr.in, species="mm10", request_interval = 300)
  return(out.great)
})
out.tb.lst <- lapply(out.great.lst, function(out.great){
  out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great), 
                                request_interval = 300)
  return(out.tb)
})

out.great.lst <- mclapply(seq(best.K), function(i){
  gr.in <- regions.range[topic.regions[[i]], ]
  out.great <- submitGreatJob(gr.in, species="mm10", request_interval = 300)
  return(out.great)
}, mc.cores = 5)
out.tb.lst <- mclapply(out.great.lst, function(out.great){
  out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great), 
                                request_interval = 300)
  return(out.tb)
}, mc.cores = 5)

