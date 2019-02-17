# Jake Yeung
# Date of Creation: 2019-02-15
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/analyze_sorted_bulk_gene_exprs.R
# Sorted bulk analyze it

rm(list=ls())

  
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(data.table)
library(tidyr)

library(topicmodels)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(umap)
library(ggrepel)
library(biomaRt)
library(igraph)  # louvain
library(Gviz)
library(GenomicRanges)

library(nnet)
library(msgl)
library(doParallel); library(foreach)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Load LDA ----------------------------------------------------------------

inmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"
meanfilt <- 10
jbin <- "TRUE"; kstr <- "15_20_25_30_35"
# jbin <- "FALSE"; kstr <- "15_20_25_30"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
jmark <- "H3K4me1"
# jmark <- "H3K9me3"
indir <- paste0("lda_outputs.meanfilt_", 
                meanfilt,
                ".cellmin_100.cellmax_500000.binarize.", jbin, 
                ".no_filt")
fname <- paste0("lda_out_meanfilt.BM-", jmark, ".CountThres0.K-", kstr, ".Robj")
inf <- file.path(inmain, indir, fname)
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

out.lda <- ChooseBestLDA(out.lda)
(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)


# Load file  --------------------------------------------------------------

# load("~/data/scchic/public_data/E-MTAB-3079-atlasExperimentSummary.Rdata", v=T)
dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Link with H3K4me1 ------------------------------------------------------------

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms


# Are top hits correlated with gene expession -----------------------------

top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

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
top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))


# Summarize in logistic regression ----------------------------------------

top.thres <- 0.999
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
top.regions <- unique(unlist(topic.regions))
top.peaks.sub <- subset(top.peaks.annotated, term %in% top.regions)
genes.keep <- unique(top.peaks.sub$SYMBOL)

print(length(genes.keep))

# subset exprs matrix and fit
dat.sub <- subset(dat.long, Gene_Name %in% genes.keep)
dat.sub$CellTypeY <- relevel(as.factor(dat.sub$CellType), ref = "fetal_liver_hematopoietic_progenitor_cell")

# jfit <- multinom(CellTypeY ~ Gene_Name * zscore, data = dat.sub, MaxNWts = 80000)

X <- model.matrix(~ Gene_Name * zscore, data = dat.sub)

Y <- dat.sub$CellTypeY
cl <- makeCluster(4)
registerDoParallel(cl)
jfit <- msgl::cv(x = X, classes = Y, standardize = FALSE, alpha = 1, lambda=0.1, use_parallel = TRUE)

saveRDS(jfit, file = paste0("~/data/scchic/robjs/multinom_fit.", top.thres, ".rds"))
