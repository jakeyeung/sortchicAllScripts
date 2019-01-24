# Jake Yeung
# Date of Creation: 2019-01-23
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_LDA_100kb_H3K4me3.R
# H3K4me3 analysis 

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
library(JFuncs)
library(umap)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

source("scripts/Rfunctions/PlotFunctions.R")

# Functions ---------------------------------------------------------------




# Load --------------------------------------------------------------------

jbin <- "TRUE"
jchip <- "H3K4me3"
# from t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/lda_analysis.h3k27me3.dropbox/lda_outputs.meanfilt_1000.merge_Windows2.cellmin_NA.cellmax_NA/lda_out.meanfilt.K-12.Robj
# from t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-5_10_15_20_25.Robj
# inf <- paste0("/private/tmp/lda_output_binned/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_10_15_20_25.Robj")
if (jbin){
  inf <- paste0("/private/tmp/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_10_15_20_25.Robj")
} else {
  inf <- paste0("/private/tmp/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_15_25.Robj")
}

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# out.lda <- out.lda[[5]]
# print(out.lda@loglikelihood)
out.lda <- ChooseBestLDA(out.lda)
(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

nn=35
nnterms <- 15
jmetric='euclidean' 
jmindist=0.15
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# check your umap settings
jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)


# Plot dat umap -----------------------------------------------------------
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))

# Plot terms umap ---------------------------------------------------------
top.thres <- 0.995
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
top.regions <- unique(unlist(topic.regions))
terms.mat <- t(tm.result$terms)[top.regions, ]
dat.umap.terms <- umap(terms.mat, config = custom.settings.terms)
# downsample rows for plotting purposes
downsamp.i <- sample(seq(nrow(dat.umap.terms$layout)), size = round(0.1 * nrow(dat.umap.terms$layout)), replace = FALSE)
jcol.rgb.terms <- lapply(seq(kchoose), ColorsByGamma, terms.mat[downsamp.i, ], c("pink", "red", "darkred"))
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb.term, jtopic){
  plot(dat.umap.terms$layout[downsamp.i, 1], dat.umap.terms$layout[downsamp.i, 2],
       col = jcol.rgb.term, pch = 20, asp = 0.75,
       main = paste("Peak Weights, T", jtopic))
}, jcol.rgb.terms, seq(kchoose))

# how many cells in the island?
cell.assign <- apply(topics.mat, 1, which.max)

top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

# annotate regions?
regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd), 
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)

hit.peaks <- subset(regions.annotated, abs(distanceToTSS) < 50000 & grepl("Hbb", SYMBOL))$region_coord

jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)


# which has Hbb?
subset(top.peaks, topic == 10)

# plot top peak
jpeak <- hit.peaks[[1]]

jpeak <- "chr6:120580000-120680000"
PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, umap.settings=umap.settings)

# find top hits 
subset(top.peaks, topic == 12)
subset(top.peaks, term %in% hit.peaks)

# subset(top.peaks, topic == 13)

subset(top.peaks, topic == 21)
subset(top.peaks, topic == 2)
subset(top.peaks, topic == 1)
subset(top.peaks, topic == 16)
subset(top.peaks, topic == 22)

# Find variable genes -----------------------------------------------------

mat.norm <- t(tm.result$topics %*% tm.result$terms)

# # find variable genes
# var.out <- GetVar(tm.result, regions.annotated)
# var.out <- var.out %>%
#   arrange(desc(Var))

# top.peaks.var <- dplyr::left_join(top.peaks, subset(var.out, select = c(region_coord, Var)), by = c("term" = "region_coord"))

jpeak <- "chr6:120580000-120680000"
jpeak <- "chr6:120580000-120680000"
jpeak <- "chr6:120580000-120680000"

PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)


# Get graph object --------------------------------------------------------




# Plot sushi plots in genomic regions -------------------------------------

# write 

