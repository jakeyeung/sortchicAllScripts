# Jake Yeung
# Date of Creation: 2019-01-22
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_analysis_H3K4me1.R
# Do downstream analysis on H3K4me1 to find clusters, interesting genes, and regulators.


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

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

PlotImputedPeaks <- function(tm.result, regions.annotated, peaks.keep, jdist, jchip, show.plot=TRUE, return.plot.only=FALSE, use.count.mat=NULL){
  # regions.sub <- subset(regions.annotated, abs(distanceToTSS) < jdist & grepl(jgene, SYMBOL)) %>%
  #   mutate(coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"),
  #          coord.kb = paste(paste(seqnames, start / 1000, sep = ":"), end / 1000, sep = "-"))
  # if (nrow(regions.sub) == 0){
  #   warning(paste("empty regions.sub for jgene", jgene, "jchip", jchip))
  #   m <- ggplot()
  #   return(m)
  # }
  # peaks.keep <- regions.sub$coord
  # peaks.keep.kb <- regions.sub$coord.kb
  peaks.keep.kb <- peaks.keep
  
  if (is.null(use.count.mat)){
    mat.norm <- t(tm.result$topics %*% tm.result$terms)
  } else {
    mat.norm <- use.count.mat
  }
  
  jlab <- peaks.keep
  
  # run on imputed matrix
  row.i <- which(rownames(mat.norm) %in% peaks.keep)
  
  if (length(row.i) > 1){
    # print(paste("Merging", length(row.i)), "rows")
    jcounts.norm <- Matrix::colSums(mat.norm[row.i, ])
  } else {
    jcounts.norm <- mat.norm[row.i, ]
  }
  # print(range(jcounts.norm))
  jcol.counts <- ColorsByCounts(jcounts.norm, nbreaks = 100, colvec = c("lightblue", "darkblue"))
  # plot topics soft clustering weights
  topics.mat <- tm.result$topics
  nn=5
  nnterms=15
  jmetric='euclidean' 
  jmindist=0.1
  custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
  custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
  dat.umap <- umap(topics.mat, config = custom.settings)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  jpeaks.str <- paste(peaks.keep.kb, collapse = ",")
  jmain <- paste0(jchip, " ", jlab, ' npeaks ', length(row.i), "\n", jpeaks.str)
  
  # prepare plot object
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2], 
                    jcol = jcol.counts,
                    counts.norm = jcounts.norm)
  m <- ggplot(dat, aes(x = umap1, y = umap2, col = jcol.counts)) + geom_point() + scale_color_identity() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    ggtitle(jmain)
  
  if (show.plot){
    print(m)
  }
  
  # plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, 
  #      main = jmain, col = jcol.counts, asp = 0.75)
  if (!return.plot.only){
    return(list('m' = m))
  } else {
    return(m)
  }
}


# Load H3K4me1 data -------------------------------------------------------

jchip <- "H3K4me1"

jdist <- 1000L
jmean <- 1
jmin <- 100L
jmax <- 500000L
# binarize <- "TRUE";  jtops <- "5_7_10_12_15_20_25_30"
binarize <- "FALSE"; jtops <- "15_20_25_30_35"

jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
infbase <- basename(inf)
infbase <- strsplit(infbase, ".Robj")[[1]][[1]]

# 0.98 threshold 
# inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
# 0.96 threshold 
inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.0.96.Robj"))

inf.mc <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jchip, ".datadir_mc_f.Rda"))

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.GREAT))
assertthat::assert_that(file.exists(inf.mc))


load(inf, v=T)
load(inf.GREAT, v=T)

# normalize count.mat by total sum
count.mat <- sweep(count.mat, MARGIN = 2, STATS = Matrix::colSums(count.mat), FUN = "/") * 10^6

lapply(topic.regions, length)


# what's the sparsity?
print(1 - Matrix::nnzero(count.mat) / length(count.mat))

regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- rownames(regions.range)

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

topics.mat <- tm.result$topics

# print(head(tm.result$topics))
# cluster the terms, but first get interesting terms
top.regions <- unique(unlist(topic.regions))  # topic.regions defined by threshold, 98 or 96th percentile of top weights in each column of the betas matrix
terms.mat <- t(tm.result$terms)[top.regions, ]


nn=5
nnterms=15
jmetric='euclidean' 
jmindist=0.1
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)


jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))


# plot terms: slow?

dat.umap.terms <- umap(terms.mat, config = custom.settings.terms)
# downsample rows for plotting purposes
downsamp.i <- sample(seq(nrow(dat.umap.terms$layout)), size = round(0.1 * nrow(dat.umap.terms$layout)), replace = FALSE)
jcol.rgb.terms <- lapply(seq(kchoose), ColorsByGamma, terms.mat[downsamp.i, ], c("pink", "darkred"))

par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb.term, jtopic){
  plot(dat.umap.terms$layout[downsamp.i, 1], dat.umap.terms$layout[downsamp.i, 2], 
       col = jcol.rgb.term, pch = 20, asp = 0.75,
       main = paste("Peak Weights, T", jtopic))
}, jcol.rgb.terms, seq(kchoose))

# Topic 17 hits? Which ones are Hbb? Can we merge peaks together? 

# these are outliers? ignore it
# outlier peaks ?? 
outlier.peaks <- names(which(dat.umap.terms$layout[, 2] < -30))
# where are they on the beta terms ?
terms.sub <- terms.mat[outlier.peaks, ]
outlier.topics <- unique(apply(terms.sub, 1, which.max))  # topic 7 
# which cells have high topic 7?
cell.assign <- apply(topics.mat, 1, which.max)

# find differential gene expression?
# use original counts matrix?

# topic 17 versus all
jtopic <- 15
count.sub <- count.mat[top.regions, ]
# weights.17 <- topics.mat[, 17]
# weights.17not <- 1 - topics.mat[, 17]
weights.in <- as.integer(cell.assign == jtopic)
weights.innot <- as.integer(cell.assign != jtopic)

# find the islands
weights.in <- as.integer(dat.umap$layout[, 2] < -3)
weights.innot <- as.integer(! weights.in)

col.binary <- weights.in + 1

# see where jtopic is in the map
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), 
     col = col.binary, asp = 0.75)


exprs.in <- apply(count.mat, 1, weighted.mean, w = weights.in)
exprs.out <- apply(count.mat, 1, weighted.mean, w = weights.innot)

# log fold change
log.fc <- log2(exprs.in + 1) - log2(exprs.out + 1)

print(head(log.fc[order(abs(log.fc), decreasing = TRUE)]))
print(head(log.fc[order(log.fc, decreasing = TRUE)]))

exprs.dat <- data.frame(exprs.in = exprs.in, exprs.out = exprs.out, log.fc = log.fc)
exprs.dat$coord <- rownames(exprs.dat)


# where is Hbb?
jsub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & grepl("Hbb", SYMBOL)) %>%
  mutate(coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"))
jpeaks <- jsub$coord

exprs.dat$jlab <- sapply(exprs.dat$coord, function(x) ifelse(x %in% jpeaks, x, ""))

# compare with top hits in beta matrix
top.hits <- sort(tm.result$terms[jtopic, ], decreasing = TRUE)

print(head(top.hits))
print(head(exprs.dat %>% arrange(desc(log.fc)), n = 50))

ggplot(exprs.dat, aes(x = log2(exprs.in + 1), y = log.fc, label = jlab)) + geom_point() + geom_text() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# plot hits on map

jpeak <- "chr9:116069000-116222000"
jpeak <- "chr17:24186000-24209000"
jpeak <- "chr5:143610000-143690000"
jpeak <- "chr12:31909000-32028000"
jpeak <- "chr10:43570000-43691000"
jpeak <- "chr5:143610000-143690000"
jpeak <- "chr17:24186000-24209000"
jpeak <- "chr6:68413000-68507000"
jpeak <- jpeaks[[1]]
jpeak <- jpeaks[[2]]
jpeak <- jpeaks[[3]]

jpeak <- "chr12:31909000-32028000"
jpeak <- "chr10:43570000-43691000"

PlotImputedPeaks(tm.result, regions.annotated, jpeak, jdist = 10000, jchip = "H3K4me1", return.plot.only = TRUE, 
                 use.count.mat = count.mat)
PlotImputedPeaks(tm.result, regions.annotated, jpeak, jdist = 10000, jchip = "H3K4me1", return.plot.only = TRUE)

# system.time(
#   plots.out <- PlotAllMarks(jgene, jchips, jdist, out.objs)
# )
