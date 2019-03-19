# Jake Yeung
# Date of Creation: 2019-03-15
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/GREAT_on_bins.R
# GREAT on bins

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

source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")

mark <- "H3K9me3"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-", mark, ".CountThres0.K-15_20_25_30.Robj")
# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-", mark, ".CountThres0.K-15_20_25_30_35.Robj")

assertthat::assert_that(file.exists(inf))

top.thres <- 0.96

load(inf, v=T)

out.lda.lst <- out.lda


Kvec <- sapply(out.lda.lst, function(x) x@k)
best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]
kchoose <- best.K
out.lda <- out.lda.lst[[which(Kvec == kchoose)]]

tmResult <- posterior(out.lda)

# names(colnames(tmResult$terms)) <- colnames(tmResult$terms)
colnames(tmResult$terms) <- gsub("chr20", "chrX", colnames(tmResult$terms))
colnames(tmResult$terms) <- gsub("chr21", "chrY", colnames(tmResult$terms))


# Run GREAT --------------------------------------------------------------

# need to assign cutoff for each peak for each topic
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
# rownames(regions.annotated) <- regions.annotated$region_coord
regions.annotated$region_coord <- names(regions.range)

i <- 2  # test
gr.in <- regions.range[names(topic.regions[[i]]), ]
out.great <- submitGreatJob(gr.in, species="mm10", request_interval = 30)
out.tb <- getEnrichmentTables(out.great, ontology=availableOntologies(out.great))

ontos <- names(out.tb)
GOdata <- out.tb[[ontos[[2]]]]

names(GOdata) <- ontos

out.tb.merge <- lapply(ontos, function(onto){
  out.tb[[onto]] %>% arrange(Hyper_Raw_PValue) %>% mutate(onto = onto)
})  %>%
  bind_rows()



onto <- "MSigDB Immunologic Signatures"
onto <- "MSigDB Perturbation"

out.tb[[onto]] %>% arrange(Hyper_Raw_PValue)


