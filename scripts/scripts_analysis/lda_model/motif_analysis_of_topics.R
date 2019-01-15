# Jake Yeung
# Date of Creation: 2019-01-14
# File: ~/projects/bs-trans/jake/scripts/analysis/motif_analysis_of_topics.R
# Downstream analysis of topics

rm(list=ls())

library(dplyr)
library(topicmodels)
library(RcisTarget)
library(feather)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

GetEnrichment <- function(peaks.sub, regions.annotated, db.inf){
  # filter out 10kb up and down tss 
  # peaks.sub <- topic.regions[[jtop]]
  regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
  geneLists <- unique(regions.sub$SYMBOL)
  motifRankings <- importRankings(db.inf)
  # 1. Calculate AUC
  motifs_AUC <- calcAUC(geneLists, motifRankings)
  # 2. Select significant motifs, add TF annotation & format as table
  motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                             motifAnnot=motifAnnotations_mgi)
  # 3. Identify significant genes for each motif
  # (i.e. genes from the gene set in the top of the ranking)
  # Note: Method 'iCisTarget' instead of 'aprox' is more accurate, but slower
  motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                     geneSets=geneLists,
                                                     rankings=motifRankings, 
                                                     nCores=1,
                                                     method="aprox")
  enrich.filt <- subset(motifEnrichmentTable_wGenes, AUC > 0.005 & rankAtMax < 20000)
  return(enrich.filt)
}



# Load a topic output -----------------------------------------------------

# jchip <- "H3K4me1"
jchips <- c("H3K27me3", "H3K4me1", "H3K4me3", "H3K9me3")

for (jchip in jchips){
  print(jchip)
  load(paste0("/private/tmp/ldaAnalysisHiddenDomains_1000/lda_outputs.meanfilt_1.cellmin_1000.cellmax_50000.binarize.FALSE/lda_out_meanfilt.PZ-BM-", jchip, ".CountThres0.K-15_20_25_30_35.Robj"), v=T)
  load(paste0("/private/tmp/ldaAnalysisHiddenDomains_1000/lda_outputs.meanfilt_1.cellmin_1000.cellmax_50000.binarize.FALSE/downstream/lda_out_meanfilt.PZ-BM-", jchip, ".CountThres0.K-15_20_25_30_35.GREAT.0.96.Robj"), v=T)
  
  db.inf <- "/Users/yeung/projects/scChiC/data/motifs/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
  
  dat <- feather(db.inf)
  
  regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                  annoDb='org.Mm.eg.db'))
  regions.annotated$region_coord <- names(regions.range)
  
  # Run Rcistarget ----------------------------------------------------------
  
  # jtop <- 1
  
  etables <- lapply(topic.regions, GetEnrichment, regions.annotated, db.inf)
  saveRDS(etables, paste0("~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/motif_enrichment_tables.", jchip, ".rds"))
  
  # lapply(etables, function(x) print(head(subset(x, select = c(motif, AUC, TF_highConf)))))
  
}
