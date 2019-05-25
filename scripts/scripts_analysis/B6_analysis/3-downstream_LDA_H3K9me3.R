# Jake Yeung
# Date of Creation: 2019-05-10
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/3-downstream_LDA_H3K9me3.R
# H3K9me3 may be different settings

rm(list=ls())

library(JFuncs)
library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)
library(ggrepel)
library(tidyr)

library(hash)
library(igraph)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load files --------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# jbin <- TRUE; kstr <- "25_30_40_50"
jbin <- FALSE; kstr <- "30_40_50"
keep.top.genes <- 150


# jmark <- "H3K4me3"
jmark <- "H3K9me3"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.", jbin, ".no_filt/lda_out_meanfilt.B6_", jmark, "_pcutoff_0.CountThres0.K-", kstr, ".Robj")
# inf <- paste09"/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_40_50.Robj"
assertthat::assert_that(file.exists(inf))

# Process LDA -------------------------------------------------------------

# kchoose <- "auto"
kchoose <- 50
out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)
# load(inf, v=T)
print(paste("K:", out.objs$out.lda@k))

# Plot UMAP ---------------------------------------------------------------

jmetric.louv='euclidean'
jmindist.louv=0.3
jseed.louv=123

# nn.louv.new <- c(150, 100, 33, 31)

jmindist.new <- c(0.5, 0.5, 0.3, 0.5)
nn.new <- c(60, 60, 45, 30)
custom.settings.new.lst <- mapply(function(x, y) GetUmapSettings(x, jmetric.louv, y, jseed.louv), nn.new, jmindist.new, SIMPLIFY = FALSE)
# custom.settings.louv.new.lst <- lapply(nn.louv.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))

names(custom.settings.new.lst) <- jmarks
# names(custom.settings.louv.new.lst) <- jmarks

dat.umap <- umap(out.objs$tm.result$topics, config = custom.settings.new.lst[[jmark]])

dat.umap.long <- data.frame(cell = rownames(dat.umap$layout), umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], mark = jmark)

# PlotXYWithColor(dat.umap.long, xvar = "umap1", yvar = "umap2")

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Do louvain clustering ---------------------------------------------------

# x <- 100
# x <- 115
x <- 60  # same louvain
jseed <- jseed.louv
jsettings <- GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed)
louv.hash <- DoLouvain(out.objs$tm.result$topics, jsettings)

# add to umap
dat.umap.long$louvain <- as.character(sapply(as.character(dat.umap.long$cell), function(x) louv.hash[[x]]))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + ggtitle(paste(x, jseed))

# After reannotating, do the clustering  ----------------------------------

tm.result <- out.objs$tm.result
# get top topics 
# analyze topic matrix across cells
topics.mat <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics))
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

terms.long <- data.frame(term = colnames(tm.result$terms), as.data.frame(t(tm.result$terms)), stringsAsFactors = FALSE) %>%
  gather(key = "topic", value = "weight", -term) %>%
  mutate(topic = gsub("X", "", topic)) %>%
  group_by(topic) %>%
  arrange(desc(weight)) %>%
  mutate(rnk = seq(length(weight))) %>%
  rowwise()

# Load facs ---------------------------------------------------------------



# Reannotate genes --------------------------------------------------------



# filter out top 1000 peaks

# annotate terms
terms.filt.top <- terms.long %>%
  filter(rnk < 1000) %>%
  rowwise()



tss.dat <- fread("/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", col.names = c("seqnames", "start", "end", "tssname"))
tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])

annots.biomart <- out.objs$regions.annotated %>% 
  mutate(midpt = start + (end - start) / 2) %>%
  filter(region_coord %in% terms.filt.top$term)


annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)

annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
  mutate(midpt = start + round(width / 2),
         midpt.1 = start.1 + round(width.1 / 2),
         dist.to.tss = midpt.1 - midpt)

# filter closest
out2.df.closest <- out2.df %>%
  group_by(region_coord) %>%
  filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

# Find interesting topics -------------------------------------------------

terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)  # error stack overflow?

terms.filt <- terms.filt.top %>%
  mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
  filter(!is.na(termgene)) %>%
  mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
  group_by(gene) %>%
  filter(weight == max(weight))




# plot --------------------------------------------------------------------


topics.all <- topics.sum$topic
names(topics.all) <- topics.all

top.genes.all.lst <- lapply(topics.all, function(jtopic){
  jtmp <- subset(terms.filt, topic == jtopic)$gene
  jtmp <- jtmp[1:min(length(jtmp), keep.top.genes)]
  return(jtmp)
})

# get averag zscore across celltypes for each topic
zscores.avg.all <- lapply(topics.all, function(jtopic){
  jtopic <- as.character(jtopic)
  top.genes <- top.genes.all.lst[[jtopic]]
  jsub <- subset(dat.long, Gene_Name %in% top.genes) %>%
    group_by(CellType) %>%
    summarise(zscore = mean(zscore)) %>%
    ungroup() %>%
    mutate(weight = exp(zscore) / sum(exp(zscore)))
  jsub$topic <- jtopic
  return(jsub)
}) %>%
  bind_rows()



# Order topics by celltype identity ---------------------------------------

zscores.sum <- zscores.avg.all %>%
  group_by(topic) %>%
  summarise(entropy = -sum(weight * log(weight))) %>%
  arrange(entropy)
  

# Plot top hits -----------------------------------------------------------

print(topics.sum)
# plot a topic
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X24")
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X13")
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X26")
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X34")

ntopicsfilt <- round(nrow(topics.sum) / 2)
# topics.keep <- topics.sum$topic[1:ntopicsfilt]
topics.keep <- zscores.sum$topic[1:ntopicsfilt]
names(topics.keep) <- topics.keep

topics.ordered <- zscores.sum$topic

genes.keep <- subset(terms.filt, topic %in% topics.keep)$gene # remove meaningless topics



top.genes.lst <- lapply(topics.keep, function(jtopic){
  jtmp <- subset(terms.filt, topic == jtopic)$gene
  jtmp <- jtmp[1:min(length(jtmp), keep.top.genes)]
  return(jtmp)
})

# do all topics



# get averag zscore across celltypes for each topic
zscores.avg <- lapply(topics.keep, function(jtopic){
  jtopic <- as.character(jtopic)
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


# Check top topics --------------------------------------------------------

# topvecs <- topics.sum$topic[1:10]

dat.withtopic <- left_join(dat.umap.long, topics.long) 

pdf(paste0("~/data/scchic/pdfs/B6_figures/B6_", jmark, "_bin_", jbin, "_k_", kchoose, "_topics_plot_", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (topvec in topics.keep){
  m1 <- PlotXYWithColor(dat.withtopic %>% filter(topic == topvec), xvar = "umap1", yvar = "umap2", cname = "weight", jtitle = topvec)
  print(m1)
}
dev.off()




# do all the topics 

# this goes into main figure probably

pdf(paste0("~/data/scchic/pdfs/B6_figures/B6_", jmark, "_bin_", jbin, "_k_", kchoose, "_celltypes_by_topic_only.keepgenes.", keep.top.genes, ".", Sys.Date(), ".winbugfix.pdf"), useDingbats = FALSE)
# ctypes.filt <- c("T_cell", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "granulocyte", "Kit_and_Sca1-positive_hematopoietic_stem_cell")
topskeep <- topics.sum$topic[1:5]
m1 <- ggplot(zscores.avg %>% filter(topic %in% topskeep), aes(y = topic, x = CellType, size = weight, color = weight)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m1)
for (jtop in topics.keep){
  # plot also the UMAP 
  m.umap <- PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = paste0("X", jtop), jtitle = jtop)
  jsub.avg <- subset(zscores.avg.all %>% filter(topic == jtop)) %>%
    arrange(desc(zscore))
  ctypevec <- jsub.avg$CellType
  jsub <- subset(dat.long, Gene_Name %in% top.genes.all.lst[[as.character(jtop)]])
  jsub.avg$CellType <- factor(jsub.avg$CellType, levels = ctypevec)
  jsub$CellType <- factor(jsub$CellType, levels = ctypevec)
  # reorder by avg?
  # sort things
  m.genes <- ggplot(jsub, aes(x = CellType, y = zscore)) + geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) + 
    ggtitle(paste("Topic", jtop))
  m.top <- ggplot(jsub.avg, aes(x = CellType, y = weight)) + geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) + 
    ggtitle(paste("Topic:", jtop))
  # show top loadings
  hits.sub <- subset(terms.filt, topic == jtop & rnk <= keep.top.genes)
  hits.sub <- OrderDecreasing(hits.sub, jfactor = "termgene", jval = "weight")
  m.hits <- ggplot(hits.sub, aes(x = termgene, y = log10(weight), label = gene)) +
    geom_text_repel() +
    theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Topic:", jtop))
  print(m.umap)
  print(m.genes)
  print(m.top)
  print(m.hits)
}
dev.off()


# Write objects -----------------------------------------------------------



