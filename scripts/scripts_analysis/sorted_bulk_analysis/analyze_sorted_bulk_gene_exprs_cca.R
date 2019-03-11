# Jake Yeung
# Date of Creation: 2019-02-21
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/analyze_sorted_bulk_gene_exprs_ttest.R
# With t-test

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

library(scales)

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
nn=40
nnterms <- 15
jmetric='euclidean' 
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# pdf(paste0("/tmp/", jmark, "_output.pdf"), useDingbats = FALSE)
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], -dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))
# dev.off()

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


# Do CCA to find cell-types that explain topics ---------------------------

beta.mat <- top.peaks.annotated %>% 
  ungroup() %>%
  mutate(topic = paste0("topic_", topic)) %>% 
  group_by(topic, SYMBOL) %>%
  filter(rnk == min(rnk)) %>%
  spread(data = ., key = "topic", value = "beta", fill = 0) %>%
  dplyr::select(-c(term, rnk))

exprs.mat <- dat.long %>%
  ungroup() %>%
  dplyr::select(c(Gene_Name, CellType, zscore)) %>%
  group_by(Gene_Name, CellType) %>%
  filter(zscore == max(zscore)) %>%
  spread(data = ., key = "CellType", value = "zscore")

# filter out top genes??

genes.keep <- intersect(unique(beta.mat$SYMBOL), unique(exprs.mat$Gene_Name))

X <- as.matrix(subset(beta.mat, select = -SYMBOL))
rownames(X) <- beta.mat$SYMBOL
X <- X[which(genes.keep %in% rownames(X)), ]
X <- t(scale(t(X), center = TRUE, scale = TRUE))
X <- scale(X, center = TRUE, scale = TRUE)

Y <- as.matrix(subset(exprs.mat, select = -Gene_Name))
rownames(Y) <- exprs.mat$Gene_Name
Y <- Y[which(genes.keep %in% rownames(Y)), ]
Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
Y <- scale(Y, center = TRUE, scale = TRUE)

out <- cancor(X, Y, xcenter = FALSE, ycenter = FALSE)


x1 <- 3
x2 <- 4

x1 <- 5
x2 <- 6

x1 <- 1
x2 <- 2


par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(out$xcoef[, x1], out$xcoef[, x2], pch = 20)
text(out$xcoef[, x1], out$xcoef[, x2], labels = rownames(out$xcoef))
plot(out$ycoef[, x1], out$ycoef[, x2], pch = 20)
text(out$ycoef[, x1], out$ycoef[, x2], labels = rownames(out$ycoef))


# Summaerize across all tpics ---------------------------------------------

topn <- 500
jtopics <- as.character(sort(unique(top.peaks.annotated$topic)))

# get average zscore 
top.genes.lst <- lapply(jtopics, function(jtopic){
  subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]
})

# get averag zscore across celltypes for each topic
zscores.avg <- lapply(as.numeric(jtopics), function(jtopic){
  top.genes <- top.genes.lst[[jtopic]]
  jsub <- subset(dat.long, Gene_Name %in% top.genes) %>%
    group_by(CellType) %>%
    summarise(zscore = mean(zscore)) %>%
    ungroup() %>%
    mutate(weight = exp(zscore) / sum(exp(zscore)))
  jsub$topic <- jtopic
  return(jsub)
}) %>%
  bind_rows()

# sort by range
zscores.ranked <- zscores.avg %>%
  group_by(topic) %>%
  summarise(diffrange = diff(range(zscore)),
            zscoremax = max(zscore)) %>%
  arrange(desc(diffrange))

zscores.ranked.bymax <- zscores.ranked %>%
  arrange(desc(zscoremax))
  
# jtopics.reordered <- zscores.ranked$topic
jtopics.reordered <- zscores.ranked.bymax$topic

# jtopics.keep <- subset(zscores.ranked, diffrange > 1.13)$topic
# jtopics.keep <- zscores.ranked$topic[1:10]
jtopics.keep <- zscores.ranked.bymax$topic[1:10]

zscores.avg$topic <- factor(zscores.avg$topic, levels = jtopics.keep)

# ggplot(zscores.avg %>% filter(topic %in% jtopics.keep), aes(y = as.character(topic), x = CellType, size = zscore, color = zscore)) + geom_point()
m.summary <- ggplot(zscores.avg %>% filter(topic %in% jtopics.keep), 
                    aes(y = topic, x = CellType, size = weight, color = weight)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

m.summary.t <- ggplot(zscores.avg %>% filter(topic %in% jtopics.keep), 
                    aes(x = topic, y = CellType, size = weight, color = weight)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)

outdir <- "~/Dropbox/scCHiC_figs/FIG4_BM/sorted_bulk_correlation"
dir.create(outdir)
# pdf(file.path(outdir, paste0(jmark, "_sorted_bulk_correlation.sortbyZmax.pdf")), useDingbats = FALSE)

dat.umap.gather <- gather(dat.umap.long, key = "topic", value = "weight", c(-umap1, -umap2, -cell))
dat.umap.gather$topic <- factor(dat.umap.gather$topic, levels = as.numeric(sort(unique(dat.umap.gather$topic))))

dat.umap.gather$colval <- log10(dat.umap.gather$weight * 10^6 + 1)
# plot all tpics
m.umap.all <- ggplot(dat.umap.gather, aes(x = umap1, y = umap2, color = colval)) + geom_point(size = 0.1) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_gradient2(low = muted('blue'), mid = "gray85", high = muted('red'), midpoint = mean(dat.umap.gather$colval)) + facet_wrap(~topic, scales = "free")
print(m.umap.all)

print(m.summary)
print(m.summary.t)


# Can we do t-test to find statistically significant cell-types? ----------

# take topic 1 as test
jtopic <- 22  # Bcells?
# jtopic <- 1
jsub <- subset(top.peaks.annotated, topic == jtopic)[1:15, ] %>%
  group_by(SYMBOL) %>%
  filter(rnk == min(rnk))

dat.sub <- subset(dat.long, Gene_Name %in% jsub$SYMBOL)
# dat.sub <- subset(dat.long, Gene_Name %in% top.genes)
dat.sub.sorted.summarised <- dat.sub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
jlevels <- as.character(dat.sub.sorted.summarised$CellType)
dat.sub$CellType <- factor(dat.sub$CellType, levels = jlevels)

dat.sub <- dat.sub %>%
  group_by(Gene_Name) %>%
  mutate(weight = exp(zscore) / sum(exp(zscore)))

genes.sub <- intersect(jsub$SYMBOL, dat.sub$Gene_Name)

# add beta values to dat.sub
dat.merge <- left_join(dat.sub, jsub, by = c("Gene_Name" = "SYMBOL")) %>%
  ungroup() %>%
  mutate(logbeta = log10(beta),
         logbeta.centered = scale(logbeta, center = TRUE, scale = TRUE),
         logbeta.weight = exp(beta) / sum(exp(beta)), 
         dummy = 1)

# fit <- lm(beta ~ zscore : CellType, data = dat.merge)
# fit <- lm(logbeta ~ zscore : CellType, data = dat.merge)
# fit <- lm(dummy ~ zscore : CellType, data = dat.merge)
# fit <- lm(beta ~ 0 + weight : CellType, data = dat.merge)
fit <- lm(logbeta.weight ~ weight : CellType, data = dat.merge)

fit.coef <- sort(coef(fit)[2:length(coef(fit))], decreasing = TRUE)
# fit.coef <- sort(coef(fit)[1:length(coef(fit))], decreasing = TRUE)

plot(fit.coef)
text(fit.coef, labels = names(fit.coef))

ggplot(dat.sub, aes(x = CellType, y = zscore)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dat.sub, aes(x = CellType, y = weight)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dat.sub %>% filter(Gene_Name == "4930515G16Rik"), aes(x = CellType, y = zscore)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

