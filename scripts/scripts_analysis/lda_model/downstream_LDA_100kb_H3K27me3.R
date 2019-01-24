# Jake Yeung
# Date of Creation: 2019-01-22
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/downstream_LDA_100kb_H3K27me3.R
# Analysis of some old analysis from H3K27me3

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



# from t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/lda_analysis.h3k27me3.dropbox/lda_outputs.meanfilt_1000.merge_Windows2.cellmin_NA.cellmax_NA/lda_out.meanfilt.K-12.Robj
inf <- "/private/tmp/lda_output_binned/lda_out.meanfilt.K-12.Robj"

load(inf, v=T)

kchoose <- out.lda@k
tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

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

# maybe topic 7 is interesting

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

hox.peaks <- subset(regions.annotated, abs(distanceToTSS) < 50000 & grepl("Hox", SYMBOL))$region_coord

# check out topic 12
subset(top.peaks, topic == 12)

subset(top.peaks, topic == 10)

subset(top.peaks, topic == 1)

subset(top.peaks, topic == 2)

# intersting topics
# Topic 12: Hox clusters
# Topic 10: Blnk, Foxo1, Pgcp1, Ebf1, Fsip1: Bcells?
# Topic 2: Prolactin region: Ezh2 Trophoblast? 

# top hits for sex chromosomes? 


# Find variable peaks -----------------------------------------------------


