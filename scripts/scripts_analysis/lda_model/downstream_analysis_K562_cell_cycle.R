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
# inf <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K27me3.CountThres0.K-10_15_20_25.Robj"

# inf <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K4me3.CountThres0.K-10_15_20_25.binarizeTrue.Robj"
# inf.GREAT <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K4me3.CountThres0.K-10_15_20_25.GREAT.Robj"

inf <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K27me3.CountThres0.K-10_15_20_25.binarizeTrue.Robj"
inf.GREAT <- "/private/tmp/K562_LDA/lda_out_meanfilt.PZ-K562-H3K27me3.CountThres0.K-10_15_20_25.GREAT.Robj"

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

phasenames <- c("G1", "G2M", "S")
phasecols <- c("red", "blue", "black")
colhash <- hash(phasenames, phasecols)
jcol.phase <- sapply(cellcycles, function(x) colhash[[x]], USE.NAMES = FALSE)

tmResult <- posterior(out.lda)

nn <- 5
# jmetric <- 'pearson'
# jmetric <- 'cosine'
jmetric <- 'euclidean'
jmindist <- 0.5
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

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col = jcol.phase)
legend(2, -4, legend=phasenames, col=phasecols, pch = 20)

dat.pca <- prcomp(tmResult$topics, center = TRUE, scale. = TRUE)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(dat.pca$x[, 1], dat.pca$x[, 2], col = jcol.phase, pch = 20)
# plot(dat.pca$x[, 2], dat.pca$x[, 3], col = jcol.phase, pch = 20)
# plot(dat.pca$x[, 3], dat.pca$x[, 4], col = jcol.phase, pch = 20)
# plot(dat.pca$x[, 4], dat.pca$x[, 5], col = jcol.phase, pch = 20)
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


# Load GREAT --------------------------------------------------------------

load(inf.GREAT, v=T)

jontos <- unique(unlist(lapply(out.tb.lst, function(x) names(x))))
cutoff <- 0.05
foldchange <- 1.5

ontology <- jontos[[2]]
order.by <- "Hyper_Fold_Enrichment"
jdecreasing <- TRUE
top <- 3

# dotplot from cisTopics
topics <- which(lapply(1:length(out.tb.lst), function(i) !is.null(out.tb.lst[[i]][[ontology]])) == TRUE)
GOdata <- lapply(topics, function(i) out.tb.lst[[i]][[ontology]])
GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by], decreasing = jdecreasing),])
# where is cell cycle
GOdata.long <- bind_rows(GOdata)
subset(GOdata.long, grepl("cycle", name)) %>% filter(Hyper_Raw_PValue < 0.05) %>% arrange(Hyper_Raw_PValue)

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

