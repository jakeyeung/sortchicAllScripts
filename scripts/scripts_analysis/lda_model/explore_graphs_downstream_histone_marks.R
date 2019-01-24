# Jake Yeung
# Date of Creation: 2019-01-16
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/explore_graphs_downstream_histone_marks.R
# Can we do better than UMAP for LDA output?

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

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")



# Parameters to get files -------------------------------------------------


jchip <- "H3K4me1"
jchip <- "H3K4me3"
jchip <- "H3K9me3"
jchip <- "H3K27me3"

jdist <- 1000L
jmean <- 1
jmin <- 100L
jmax <- 500000L
binarize <- "TRUE";  jtops <- "15_20_25_30_35"
# binarize <- "FALSE"; jtops <- "15_20_25_30_35"

jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
infbase <- basename(inf)
infbase <- strsplit(infbase, ".Robj")[[1]][[1]]

# 0.98 threshold 
# inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
# 0.96 threshold 
inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.0.96.Robj"))

load(inf, v=T)
load(inf.GREAT, v=T)

# Get best K --------------------------------------------------------------

lapply(topic.regions, length)


# what's the sparsity?
print(1 - Matrix::nnzero(count.mat) / length(count.mat))

if (exists("regions.range")){
  regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                  annoDb='org.Mm.eg.db'))
  regions.annotated$region_coord <- rownames(regions.range)
}


# if already loaded the GREAT object, then don't need to choose best K

if (length(out.lda) > 1){
  out.lda.lst <- out.lda
  # we did multicore, so which one to choose?
  Kvec <- sapply(out.lda.lst, function(x) x@k)
  best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]
  # plot loglikelihood
  par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood), 'o')
  kchoose <- best.K
  out.lda <- out.lda.lst[[which(Kvec == kchoose)]]
  print(paste("Likelihood: ", out.lda@loglikelihood))
} else {
  kchoose <- out.lda@k
}

tm.result <- posterior(out.lda)

mat.norm <- t(tm.result$topics %*% tm.result$terms)


# Get variable peaks ------------------------------------------------------

# find peaks with largest range
# maxmin <- sort(apply(mat.norm, 1, function(jcol) abs(diff(c(quantile(jcol, 0.8), quantile(jcol, 0.2))))), decreasing = TRUE)
maxmin <- sort(apply(mat.norm, 1, function(jcol) mad(jcol)), decreasing = TRUE)

# what are the sizes?
peak.size <- sapply(names(maxmin), function(x) ParseCoord(x)$end - ParseCoord(x)$start)

plot(seq(length(peak.size)), peak.size, pch = 20, main='MAD')

# normalize maxmin by peaksize
maxmin.norm <- maxmin / peak.size

plot(peak.size, maxmin, pch = 20)

plot(peak.size, maxmin.norm, pch = 20)

dat.var <- data.frame(Var = maxmin, peak.size = peak.size, coord = names(maxmin))

ggplot(dat.var, aes(x = peak.size, y = maxmin)) + geom_point(alpha = 0.2) 

print(head(subset(dat.var, peak.size < 10000) %>% arrange(desc(Var))))

regions.annotated <- dplyr::left_join(regions.annotated, dat.var, by = c("region_coord"="coord"))


# Plot variable regions ---------------------------------------------------



# Plot LDA map, color by Hbb region in space ------------------------------

# jgene <- "Cebpe"
jgene <- "Hox"
jgene <- "Hbb"
jgene <- "Car1"
jgene <- "Klf14"
regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 100000 & grepl(jgene, SYMBOL)) %>%
  mutate(coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"))
peaks.keep <- regions.sub$coord
jlab <- jgene

# run on imputed matrix
row.i <- which(rownames(mat.norm) %in% peaks.keep)

if (length(row.i) > 1){
  jcounts.norm <- Matrix::colSums(mat.norm[row.i, ])
} else {
  jcounts.norm <- mat.norm[row.i, ]
}


print(range(jcounts.norm))

jcol.counts <- ColorsByCounts(jcounts.norm, nbreaks = 100, colvec = c("lightblue", "darkblue"))

# plot topics soft clustering weights
topics.mat <- tm.result$topics
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)

nn=5
nnterms=15
jmetric='euclidean' 
jmindist=0.1
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
dat.umap <- umap(topics.mat, config = custom.settings)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, 
     main = paste("Normalized counts", jlab), col = jcol.counts, asp = 0.75)
print(jlab)


# Plot for all 4 marks ----------------------------------------------------



# Compare with just normalized UMI counts ---------------------------------



# par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
# mapply(function(jcol.rgb, jtopic){
#   plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
# }, jcol.rgbs, seq(kchoose))

# plot topics, plot by Hox counts (size normalized)
# Color by normalized UMI counts

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

# Find clusters and gradients ---------------------------------------------

# filter for top regions, and transpose
top.thres <- 0.99
topic.regions <- lapply(seq(out.lda@k), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})

topics.regions.common <- unique(unlist(topic.regions))
print(length(topics.regions.common))

mat.norm.sub <- mat.norm[topics.regions.common, ]

# renormalized
mat.norm.sub <- scale(mat.norm.sub, center = TRUE, scale = TRUE)

out <- umap(mat.norm.sub)

# plot(out$layout[, 1], out$layout[, 2])

# save matrix to load in python??
write.table(t(mat.norm.sub), file = paste0("outputs_R/lda_output/chip.", jchip, ".normalized_counts.txt"), quote = FALSE, sep = "\t")

