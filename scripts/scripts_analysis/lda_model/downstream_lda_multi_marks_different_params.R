# Jake Yeung
# Date of Creation: 2019-01-07
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_lda_multi_marks_different_params.R
# Run LDA on different parameters

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

pcut <- 0.5
jdist <- 1000L
jmean <- 1
jmin <- 100L
jmax <- 500000L
jchip <- "H3K4me3"
jchip <- "H3K27me3"
jchip <- "H3K9me3"
jchip <- "H3K4me1"

jdir <- paste0('/tmp/ldaAnalysisBroadpeaks_', pcut, '_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax)
inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-5_7_10_12_15_20_25_30.Robj'))
infbase <- basename(inf)
infbase <- strsplit(infbase, ".Robj")[[1]][[1]]
inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.GREAT))
load(inf, v=T)

# what's the sparsity?
print(1 - Matrix::nnzero(count.mat) / length(count.mat))

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
par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood), 'o')

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
jcounts <- Matrix::colSums(count.mat)
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

(x <- terms(out.lda, 6000))
 
# Hox genes?
hoxa <- "chr16:52"
hoxb <- "chr11:96"
hoxc <- "chr15:102"
hoxc2 <- "chr15:103"
hoxd <- "chr2:74"
hbb <- "chr7:103"

greplst <- list(hoxa, hoxb, hoxc, hoxc2, hoxd, hbb)
x.match <- lapply(greplst, function(grepstr) grep(grepstr, x, value = TRUE))
print(x.match)

# Hbb locus?


# Look at GREAT output ----------------------------------------------------

load(inf.GREAT, v=T)

jtop <- 1
jontos <- unique(unlist(lapply(out.tb.lst, function(x) names(x))))
cutoff <- 0.05
foldchange <- 1.5

ontology <- jontos[[14]]
order.by <- "Hyper_Fold_Enrichment"
jdecreasing <- TRUE
top <- 5

# dotplot from cisTopics
topics <- which(lapply(1:length(out.tb.lst), function(i) !is.null(out.tb.lst[[i]][[ontology]])) == TRUE)
GOdata <- lapply(topics, function(i) out.tb.lst[[i]][[ontology]])
GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by], decreasing = jdecreasing),])
GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][1:top,])

# annotate each GOdata by topic number
for (i in topics){
  GOdata[[i]]$topic <- i
}

GOdata <- dplyr::bind_rows(GOdata) %>%
  mutate(topic = factor(topic, levels = sort(unique(topic))))

ggplot(GOdata, aes(x = topic, y = name, size = Hyper_Fold_Enrichment, color = -log10(Hyper_Adjp_BH))) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.minor = element_blank(), axis.text=element_text(size=6)) + 
  xlab("") + ylab("")

