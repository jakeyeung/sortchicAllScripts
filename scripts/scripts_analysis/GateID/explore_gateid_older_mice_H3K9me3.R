# Jake Yeung
# Date of Creation: 2019-04-13
# File: ~/projects/scchic/scripts/scripts_analysis/GateID/explore_gateid_older_mice.R
# Explore LDA output 

rm(list=ls())

library(dplyr)
library(ggplot2)
library(topicmodels)
library(tidytext)
library(JFuncs)
library(umap)
library(tidyr)
library(data.table)
library(ggrepel)
library(hash)
library(igraph)

library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(biomaRt)


source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")

jmark <- "H3K4me9"

# Init paths and check ----------------------------------------------------


inf <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_GateID/LDA_outputs_all_GateID/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.", jmark, "_GateID.CountThres0.K-15_20_25_30_35.Robj")
assertthat::assert_that(file.exists(inf))
inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"

# Load data ---------------------------------------------------------------

# inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_GateID/lda_out_meanfilt.H3K4me1_GateID_LDA45.Robj"


pdfname <- paste0(basename(dirname(inf)), ".pdf")
# pdfdir <- "/Users/yeung/data/scchic/pdfs"
pdfdir <- "/tmp"
pdfout <- file.path(pdfdir, pdfname)


load(inf, v=T)

out.lda <- out.lda[[1]]

out.lda@k


# Load bulk sorted data ---------------------------------------------------


dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))


# Plot output -------------------------------------------------------------

tm.result <- posterior(out.lda)

# do UMAP first

custom.settings <- umap.defaults
custom.settings$n_neighbors <- 15
custom.settings$min_dist <- 0.4
custom.settings$metric <- "euclidean"


umap.out <- umap(tm.result$topics, config = custom.settings)

umap.out.long <- data.frame(umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], cell = rownames(umap.out$layout))

# Make plots --------------------------------------------------------------

pdf(pdfout, useDingbats = FALSE)


# ggplot(umap.out.long, aes(x = umap1, y = umap2)) + geom_point() + 
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

betas.long <- as.data.frame(tm.result$topics) %>% 
  mutate(cell = as.character(rownames(tm.result$topics))) %>%
  tidyr::gather(key = topic, value = beta, -cell)

# add umap to this
umap.out.long <- left_join(betas.long, umap.out.long)

m.out <- ggplot(umap.out.long, aes(x = -1 * umap1, y = -1 * umap2, color = beta)) + geom_point() + facet_wrap(~topic) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.out)


# Do louvain  -------------------------------------------------------------


jmetric.louv='euclidean'
jmindist.louv=0.2
# nn.louv <- 28
nn.louv <- 50
jseed.louv=123

custom.settings.louv <- GetUmapSettings(nn = nn.louv, jmetric = jmetric.louv, seed = jseed.louv, jmindist = jmindist.louv)

clstr.hash <- DoLouvain(tm.result$topics, custom.settings.louv, dat.umap.long = NULL)

umap.out.long$louvain <- sapply(umap.out.long$cell, function(x) clstr.hash[[x]])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")
m1 <- ggplot(umap.out.long, aes(x = -1*umap1, y = -1*umap2, color = as.character(louvain))) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
print(m1)
# save to output 
saveRDS(umap.out.long %>% mutate(umap1 = -1 * umap1, umap2 = -1 * umap2), file = "~/data/scchic/robjs/GateID_scChICseq_LouvainClusters.rds")

# Sort by most relevant topics  -------------------------------------------

# analyze topic matrix across cells
topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>%
  gather(key = "topic", value = "weight", -cell) %>%
  rowwise() %>%
  mutate(topic = as.numeric(substr(topic, 2, nchar(topic)))) %>%
  group_by(topic) %>%
  mutate(zscore = scale(weight, center = TRUE, scale = TRUE))

topics.sum <- topics.long %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(zscore < quantile(zscore, 0.97)) %>%
  mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)
print(topics.sum)

jtopics <- topics.sum$topic

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
barplot(height = topics.sum$entropy, names.arg = topics.sum$topic, xlab = "Topic", ylab = "Entropy Measure H", main = "Topics are now ordered by increasing entropy.", las = 2)






# 
# # Just do PCA -------------------------------------------------------------
# 
# pca.out <- prcomp(tm.result$topics, center = TRUE, scale. = TRUE)
# 
# plot(pca.out$x[, 1], pca.out$x[, 2], pch = 20)
# plot(pca.out$x[, 2], pca.out$x[, 3], pch = 20)
# plot(pca.out$x[, 3], pca.out$x[, 4], pch = 20)




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

top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))


# Layer on bulk expression data -------------------------------------------

topn <- 150
topn.plot <- 50
for (jtopic in jtopics){
  top.genes <- subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]
  jsub <- subset(dat.long, Gene_Name %in% top.genes)
  jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
  jlevels <- as.character(jsub.sorted.summarised$CellType)
  jsub$CellType <- factor(jsub$CellType, levels = jlevels)
  
  
  m.out <- ggplot(umap.out.long %>% filter(topic == jtopic), aes(x = -1 * umap1, y = -1 * umap2, color = beta)) + geom_point() + facet_wrap(~topic) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.out)
  
  m.ctype <- ggplot(jsub, 
                    aes(x = CellType , y = zscore)) + 
    geom_boxplot() +
    # geom_violin() +
    geom_jitter(width = 0.1, size = 0.5) +
    # geom_line() + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(jtopic)
  print(m.ctype)
  

  
  # plot top bins?
  top.sub <- top.peaks.annotated %>% filter(topic == jtopic & rnk < topn.plot) %>% arrange(rnk)
  top.sub <- OrderDecreasing(top.sub, jfactor = "term", jval = "beta")
  m.tops <- ggplot(top.sub, aes(x = term, y = beta, label = SYMBOL)) + 
    geom_text_repel() + 
    theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(paste("Topic:", jtopic))
  print(m.tops)
}

dev.off()

# check
jtest <- readRDS("/Users/yeung/data/scchic/robjs/GateID_scChICseq_LouvainClusters.rds")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")

# separate Live vs Gate
jtest$experi <- sapply(jtest$cell, function(x) substr(strsplit(x, "_")[[1]][[3]], start=1, stop = 4))
ggplot(jtest, aes(x = umap1, y = umap2, color = as.character(experi))) + geom_point(alpha = 0.25) + scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jtest, aes(x = umap1, y = umap2, color = as.character(experi))) + geom_point() + scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~experi)
