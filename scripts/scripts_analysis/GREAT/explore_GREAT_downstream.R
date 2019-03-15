# Jake Yeung
# Date of Creation: 2019-03-15
# File: ~/projects/scchic/scripts/scripts_analysis/GREAT/explore_GREAT_downstream.R
# Explore GREAT downstream

# load GTEAT Output 

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
library(JFuncs)
library(cowplot)

source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")

# Load stuff --------------------------------------------------------------

mark <- "H3K27me3"
mark <- "H3K9me3"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.TRUE.no_filt/downstream/lda_out_meanfilt.BM-", 
              mark, ".CountThres0.K-15_20_25_30.GREAT.0.96.Robj")

assertthat::assert_that(file.exists(inf))

load(inf, v=T)


# Get ontos ---------------------------------------------------------------

ontos <- names(out.tb.lst[[1]])
topics <- seq(length(out.tb.lst))

# merge everything together
out.merge <- lapply(topics, function(topic){
  out.tb <- out.tb.lst[[topic]]
  jtmp <- lapply(ontos, function(onto){
    jtmp <- out.tb[[onto]]
    jtmp$topic <- topic
    jtmp$onto <- onto
    return(jtmp)
  }) %>% bind_rows()
  return(jtmp)
}) %>% bind_rows()

jtop <- 2
jtop <- 6
jtop <- 7
jtop <- 22
jtop <- 25
jtop <- 30

lapply(ontos, function(jonto) head(out.merge %>% filter(topic == jtop & onto == jonto) %>% dplyr::select(name, Binom_Adjp_BH, topic, onto), n = 5))


# summarize MSig pathway across all cells 
print(ontos)
jonto <- "GO Cellular Component"
jonto <- "PANTHER Pathway"
jonto <- "GO Biological Process"
jonto <- "MSigDB Perturbation"
jonto <- "MSigDB Pathway"
m <- ggplot(out.merge %>% filter(onto == jonto) %>% group_by(topic) %>% top_n(n = -3, wt = Binom_Adjp_BH), 
       aes(x = as.character(topic), y = name, color = -log10(Binom_Adjp_BH), size = -log10(Binom_Adjp_BH))) + geom_point()  +
  theme_bw() + 
  ggtitle(jonto)
print(m)

# which regions are associated with Wnt signaling pathway?
jname <- "Wnt signaling pathway"
subset(out.merge, name == jname & onto == "PANTHER Pathway") %>% arrange(Hyper_Adjp_BH)
