# Jake Yeung
# Date of Creation: 2019-01-04
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/lda_binarized_topicmodels.R
# Run LDA on topic models using topic models package

rm(list=ls())

library(topicmodels)

source("scripts/Rfunctions/MetricsLDA.R")

# Load data ---------------------------------------------------------------

load("/private/tmp/lda_output/PZ-BM-H3K27me3.merged.NoCountThres.Robj", v=T)

rownames(out.lda@gamma) <- out.lda@documents
cells.keep <- out.lda@documents
peaks.keep <- out.lda@terms
cells.keep.i <- which(colnames(count.mat) %in% cells.keep)
peaks.keep.i <- which(rownames(count.mat) %in% peaks.keep)
count.mat <- count.mat[peaks.keep.i, cells.keep.i]
# remove chromoM?
count.mat <- count.dat$counts
M.peaks <- grep("chrM", rownames(count.mat), value=TRUE)
count.mat <- count.mat[which(!rownames(count.mat) %in% M.peaks), ]

load("/private/tmp/lda_outputs.meanfilt_1.merge_1000_NoM_binarize.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-20.Robj", v=T)

out.lda.lst <- out.lda

# we did multicore, so which one to choose?
Kvec <- sapply(out.lda.lst, function(x) x@k)

# plot best fit?
result <- CalculateMetrics(out.lda.lst, Kvec)

print(result)

FindTopicsNumber_plot(result)

kchoose <- 10
out.lda <- out.lda.lst[[which(Kvec == kchoose)]]

# Analyze gamma matrix ----------------------------------------------------

tmResult <- posterior(out.lda)

nn <- 10
jmetric <- 'cosine'
jmindist <- 0.001
custom.settings <- umap.defaults
custom.settings$n_neighbors <- nn
custom.settings$metric <- jmetric
custom.settings$min_dist <- jmindist

# dat.umap <- umap(out.lda@gamma, config = custom.settings)
# rownames(dat.umap$layout) <- rownames(out.lda@gamma)
dat.umap <- umap(tmResult$topics, config = custom.settings)
rownames(dat.umap$layout) <- rownames(tmResult$topics)

jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# par(mfrow=c(1,1), mar=c(1,1,1,1))
# plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's')


# color by loadings on Kvec
ColorsByGamma <- function(topic){
  # jcol <- out.lda@gamma[, topic]
  jcol <- tmResult$topics[, topic]
  colorPal <- grDevices::colorRampPalette(c("pink", "red", "darkred"))
  jcol.rgb <- colorPal(200)[as.numeric(cut(jcol,breaks = 200))]
  return(jcol.rgb)
}

jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma)

par(mfrow=c(2,5), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 1)
}, jcol.rgbs, seq(kchoose))


x <- terms(out.lda, 500)

hoxa <- "chr16:52"
hoxb <- "chr11:96"
hoxc <- "chr15:102"
hoxc2 <- "chr15:103"
hoxd <- "chr2:74"

greplst <- list(hoxa, hoxb, hoxc, hoxc2, hoxd)

x.match <- lapply(greplst, function(grepstr) grep(grepstr, x, value = TRUE))


