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

# plot example
jgene <- "Sox6"
jgene <- "Inpp4b"
jgene <- "Fam46a"
jgene <- "Hbb-y"
jgene <- "Hbb-bt"
ggplot(subset(dat.long, Gene_Name == jgene), aes(x = CellType , y = FPKM)) + 
  geom_point() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jgene)


# Link with H3K4me1 ------------------------------------------------------------

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

# settings for UMAP
nn=30
nnterms <- 15
jmetric='euclidean' 
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

pdf(paste0("/tmp/", jmark, "_output.pdf"), useDingbats = FALSE)
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], -dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))
dev.off()

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

# topic 1:
# jtopic <- 21
# jtopic <- 21
# jtopic <- 22
jtopic <- 13
jtopic <- 22
jtopic <- 26
jtopic <- 10
jtopic <- 21
jtopic <- 1

jtopic <- 10
jtopic <- 23

jtopic <- 22
jtopic <- 4
jtopic <- 30
jtopic <- 7
jtopic <- 29
jtopic <- 6
jtopic <- 26
jtopic <- 2
jtopic <- 22
jtopic <- 5

jtopic <- 1

jtopic <- 22
jtopic <- 27
jtopic <- 13

topn <- 500
# top.genes.freq <- sort(table(subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]), decreasing=TRUE)
# top.genes <- names(top.genes.freq)
top.genes <- subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]
jsub <- subset(dat.long, Gene_Name %in% top.genes)
jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
jlevels <- as.character(jsub.sorted.summarised$CellType)
jsub$CellType <- factor(jsub$CellType, levels = jlevels)
ggplot(jsub, 
       aes(x = CellType , y = zscore)) + 
  geom_boxplot() +
  # geom_violin() +
  geom_jitter(width = 0.1, size = 0.5) +
  # geom_line() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(jtopic)
  
jgene <- "Gypa"
jgene <- "Sox6"
jgene <- "Erdr1"
jgene <- "Mid1"
jgene <- "Erdr1"
ggplot(subset(dat.long, Gene_Name == jgene), aes(x = CellType , y = logFPKM)) +
  geom_point() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(jgene)
ggplot(subset(dat.long, Gene_Name == jgene), aes(x = CellType , y = zscore)) +
  geom_point() + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(jgene)

jsub <- subset(dat.long, Gene_Name == jgene)
jzscore <- max(jsub$zscore)
jsub2 <- subset(dat.long, Gene_Name == jgene) %>% 
  mutate(probC = zscore / sum(zscore),
         probCexp = exp(zscore) / sum(exp(zscore)))

dat.long %>% group_by(CellType) %>% filter(Gene_Name %in% top.genes) %>% ggplot(aes(x = zscore)) + geom_histogram() + facet_wrap(~CellType) + theme_classic() + geom_vline(xintercept = 2.25)


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
 
library(nnet)
library(msgl)
jfit <- multinom(CellTypeY ~ Gene_Name * zscore, data = dat.sub, MaxNWts = 80000)

X <- model.matrix(~ Gene_Name * zscore, data = dat.sub)

Y <- dat.sub$CellTypeY
library(doParallel); library(foreach)
cl <- makeCluster(4)
registerDoParallel(cl)
jfit <- msgl::cv(x = X, classes = Y, standardize = FALSE, alpha = 1, lambda=0.1, use_parallel = TRUE)
# test <- multinom(prog2 ~ ses + write, data = ml)
# jfit <- glm(CellTypeY ~ Gene_Name * zscore, data = dat.sub, family = multinomial)

# predict 
qplot(as.numeric(dat.sub$CellTypeY), as.numeric(predict(jfit))) + geom_jitter(width = 0.5, height = 0.5) + theme_bw() 

# accuracy
length(which(dat.sub$CellTypeY == predict(jfit))) / length(dat.sub$CellTypeY)

# Logistic regression? ----------------------------------------------------
# 
# jsub.celltype <- jsub
# jfit <- glm(CellType ~ Gene_Name * zscore, data = jsub.celltype, family = binomial())
# predict(jfit)
# 
# jtopic <- 1
# betas.sub <- subset(top.peaks.annotated, topic == jtopic) %>%
#   mutate(logbeta = log10(beta),
#          logbeta.norm = (logbeta - min(logbeta)) / (max(logbeta) - min(logbeta))) %>%
#   group_by(SYMBOL) %>%
#   filter(logbeta.norm == max(logbeta.norm))
# # predict with new values
# pred.dat <- 
# 
