# Jake Yeung
# Date of Creation: 2019-04-30
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/celltype_analysis_pretty.R
# Celltype analysis pretty
# 


rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(umap)
library(tidytext)
library(hash)

library(flexmix)


# library(devtools)
# dev_mode(T)
# install_local("/Users/yeung/projects/cellassign", force=TRUE, upgrade = "never")
# library(cellassign)


library(nnet)
library(msgl)
library(doParallel); library(foreach)

library(GenomicRanges)

library(biomaRt)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")
source("scripts/Rfunctions/Aux.R")


jmark <- "H3K4me1"

jfac <- 10^6
jpseudo <- 0

# Get trajs mixed ---------------------------------------------------------

trajs.mixed.out <- GetTrajMixed()
trajs.mixed <- trajs.mixed.out$trajs.mixed
dat.umap.mixed <- trajs.mixed.out$dat.umap.mixed


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load UMAP ---------------------------------------------------------------

inf.umap <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.umap))
load(inf.umap, v=T)

# Load exprs  -------------------------------------------------------------

inf.rdata <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
load(inf.rdata, v=T)


# Load TSS genes ----------------------------------------------------------


tssdist <- 50000
jdate <- "2019-04-22"
Kvec <- "50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, "_GeneTSS.Dedup.", jdate, ".", tssdist, ".Robj")
# inf.lda <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-50_GeneTSS.Dedup.2019-04-21.20000.Robj"

assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

out.lda <- out.lda[[length(out.lda)]]
out.lda@documents <- SwitchColnames(unname(out.lda@documents), jsplit = "-")

tss.dat <- data.frame(region_coord = sapply(out.lda@terms, function(x) strsplit(x, ";")[[1]][[1]]),
                      gene = sapply(out.lda@terms, function(x) strsplit(x, ";")[[1]][[2]]), 
                      stringsAsFactors = FALSE)
tss.dat$seqnames <- sapply(tss.dat$region_coord, GetChromo)
tss.dat$start <- sapply(tss.dat$region_coord, GetStart)
tss.dat$end <- sapply(tss.dat$region_coord, GetEnd)


# Get locations and send to biomart ---------------------------------------

jmark <- "H3K4me1"
annots.biomart <- annots.lst[[jmark]] %>% 
  mutate(midpt = start + (end - start) / 2)

annots.gr <- makeGRangesFromDataFrame(annots.lst[[jmark]] %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)

annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
  mutate(midpt = start + round(width / 2),
         midpt.1 = start.1 + round(width.1 / 2),
         dist.to.tss = midpt.1 - midpt)



# Select closest ----------------------------------------------------------

out2.df.closest <- out2.df %>%
  group_by(region_coord) %>%
  filter(dist.to.tss == min(abs(dist.to.tss))) %>%
  group_by(gene) %>%
  filter(dist.to.tss == min(abs(dist.to.tss)))

terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)


# After reannotating, do the clustering  ----------------------------------

tm.result <- tm.result.lst[[jmark]]
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

# filter out top 1000 peaks

# annotate terms
terms.filt <- terms.long %>%
  filter(rnk < 1000) %>%
  rowwise() %>%
  mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
  filter(!is.na(termgene)) %>%
  mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
  group_by(gene) %>%
  filter(weight == max(weight))


# plot a topic
# topic 22 is NK cells? 
PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = "X22")

# which topics shoul we keep?? 

ggplot(topics.sum %>% mutate(topic = factor(topic, levels = topic)), aes(x = topic, y = entropy)) + geom_point()  # note topic 22 is NK cells


ntopicsfilt <- 16
topics.keep <- topics.sum$topic[1:ntopicsfilt]
names(topics.keep) <- topics.keep

# genes.keep <- sapply(terms.filt$termgene, function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)

# genes.keep <- terms.filt$gene  # keep all
genes.keep <- subset(terms.filt, topic %in% topics.keep)$gene # remove meaningless topics



# Check top topics --------------------------------------------------------

topvecs <- topics.sum$topic[1:10]

dat.withtopic <- left_join(dat.umap.mixed[[1]], topics.long) 

for (topvec in topvecs){
  m1 <- PlotXYWithColor(dat.withtopic %>% filter(topic == topvec), xvar = "umap1", yvar = "umap2", cname = "weight", jtitle = topvec)
  print(m1)
}

# take top N genes and cluster 
keepn <- 50  # for now


# Do some celltyping by topic ---------------------------------------------


# summarize the expression of each center with a value from zscore of gene expression

keep.top.genes <- 150

top.genes.lst <- lapply(topics.keep, function(jtopic){
  jtmp <- subset(terms.filt, topic == jtopic)$gene
  jtmp <- jtmp[1:min(length(jtmp), keep.top.genes)]
  return(jtmp)
})

# do all topics

topics.all <- topics.sum$topic
names(topics.all) <- topics.all
top.genes.all.lst <- lapply(topics.all, function(jtopic){
  jtmp <- subset(terms.filt, topic == jtopic)$gene
  jtmp <- jtmp[1:min(length(jtmp), keep.top.genes)]
  return(jtmp)
})

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

# get averag zscore across celltypes for each topic
zscores.avg.all <- lapply(topics.keep, function(jtopic){
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

# do all the topics 

# this goes into main figure probably

pdf(paste0("~/data/scchic/pdfs/celltyping_analysis/celltypes_by_topic_only.keepgenes.", keep.top.genes, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
# ctypes.filt <- c("T_cell", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "granulocyte", "Kit_and_Sca1-positive_hematopoietic_stem_cell")
topskeep <- topics.sum$topic[1:5]
m1 <- ggplot(zscores.avg %>% filter(topic %in% topskeep), aes(y = topic, x = CellType, size = weight, color = weight)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m1)
for (jtop in topics.keep){
  # plot also the UMAP 
  m.umap <- PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = paste0("X", jtop), jtitle = jtop)
  jsub.avg <- subset(zscores.avg %>% filter(topic == jtop)) %>%
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
  print(m.genes)
  print(m.top)
  print(m.umap)
}
dev.off()




