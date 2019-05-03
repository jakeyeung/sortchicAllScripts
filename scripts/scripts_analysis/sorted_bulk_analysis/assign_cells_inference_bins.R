# Jake Yeung
# Date of Creation: 2019-04-26
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/assign_cells_inference_bins.R
# Bins


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


# Find marker genes -------------------------------------------------------

markers <- dat.long %>%
  mutate(is.marker = ifelse(zscore > 1.3, TRUE, FALSE))

# take markers
markers.keep <- markers %>%
  group_by(Gene_Name) %>%
  filter(max(logFPKM) > 5) %>%
  filter(any(is.marker))

print(length(unique(markers.keep$Gene_Name)))
nctypes <- length(unique(dat.long$CellType))


# Find Gaussians ----------------------------------------------------------
# rename terms

# make gene exprs matrix
dat.impute <- tm.result$topics %*% tm.result$terms
cnames.old <- unname(colnames(dat.impute))
cnames.new <- sapply(cnames.old, function(x){
  xnew <- terms.hash[[x]]
  ifelse(is.null(xnew), NA, xnew)
}) 
cnames.i <- which(!is.na(cnames.new))
dat.impute.sub <- dat.impute[, cnames.i]

colnames(dat.impute.sub) <- sapply(colnames(dat.impute.sub), function(x) terms.hash[[x]])
assertthat::assert_that(all(!is.null(colnames(dat.impute.sub))))


# Filter genes ------------------------------------------------------------

# must be in marker genes and in top 1000 of peaks?
markers <- unique(markers.keep$Gene_Name)

genes.filt2 <- genes.keep[which(genes.keep %in% markers)]  # ~700 if use all topics, ~550 if use meaningful topics, means we removed 200 noisy genes?
print(length(genes.filt2))

termgene.filt2 <- unlist(terms.filt %>%
  ungroup() %>%
  filter(gene %in% genes.filt2) %>%
  dplyr::select(termgene), use.names = FALSE)

# get top genes

dat.impute.sub2 <- dat.impute.sub[, termgene.filt2]

# do clustering by flexmix
inmat <- as.matrix(dat.impute.sub2)
# log transform the inmat and then center rows

inmat.norm <- log(inmat * 10^6)
inmat.norm <- sweep(inmat.norm, MARGIN = 1, STATS = rowMeans(inmat.norm), FUN = "-")

system.time(
  mixout <- flexmix(inmat.norm ~ 1, k = nctypes, model = FLXMCmvnorm(), control = list(classify = "weighted"))
)
mixout.post <- data.frame(cell = rownames(inmat.norm), as.data.frame(posterior(mixout)))

# plot the clustering outputs
dat.clstring <- left_join(dat.umap.mixed[[jmark]], mixout.post)

# a good one 

pdf(paste0("~/data/scchic/pdfs/celltyping_analysis/celltypes_flexmix.topicsfilt.", ntopicsfilt, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (i in seq(mixout@k)){
  jcname <- paste("V", i, sep = "")
  print(jcname)
  m1 <- PlotXYWithColor(dat.clstring, xvar = "umap1", yvar = "umap2", cname = jcname, jtitle = jcname)
  print(m1)
}
dev.off()

# summarize the expression of each center with a value from zscore of gene expression

keep.top.genes <- 250

top.genes.lst <- lapply(topics.keep, function(jtopic){
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

# do all the topics 

# this goes into main figure probably

pdf(paste0("~/data/scchic/pdfs/celltyping_analysis/celltypes_by_topic_only.keepgenes.", keep.top.genes, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
# ctypes.filt <- c("T_cell", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "granulocyte", "Kit_and_Sca1-positive_hematopoietic_stem_cell")
topskeep <- topics.sum$topic[1:5]
m1 <- ggplot(zscores.avg %>% filter(topic %in% topskeep), aes(y = topic, x = CellType, size = weight, color = weight)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m1)
for (jtop in topics.sum$topic){
  # plot also the UMAP 
  m.umap <- PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = paste0("X", jtop), jtitle = jtop)
  m.top <- ggplot(zscores.avg %>% filter(topic == jtop), aes(x = CellType, y = weight)) + geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) + 
    ggtitle(jtop)
  print(m.top)
  print(m.umap)
}
dev.off()


ggplot(zscores.avg %>% filter(topic == 26), aes(x = CellType, y = weight)) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))
ggplot(zscores.avg %>% filter(topic == 24), aes(x = CellType, y = zscore)) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))

ggplot(zscores.avg %>% filter(topic == 22), aes(x = CellType, y = weight)) + geom_bar(stat = "identity")
ggplot(zscores.avg %>% filter(topic == 22), aes(x = CellType, y = zscore)) + geom_bar(stat = "identity")
ggplot(zscores.avg %>% filter(CellType == "T_cell"), aes(x = topic, y = zscore)) + geom_bar(stat = "identity")


PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = "X26")  # erythryocytes
PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = "X22")
PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = "X24")  # BCell Lymphocytes
PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = "X10")
PlotXYWithColor(left_join(dat.umap.mixed[[jmark]], topics.mat), xvar = "umap1", yvar = "umap2", cname = "X23")



# Cool now do it on the centers -------------------------------------------

centers <- flexmix::parameters(mixout)
centers <- centers[grepl("^center", rownames(centers)), ]
centers.long <- data.frame(termgene = rownames(centers), as.data.frame(centers), stringsAsFactors = FALSE)
centers.long$gene <- sapply(centers.long$termgene, function(x) strsplit(x, ";")[[1]][[2]])


center.weights <- centers.long$Comp.3  # from PDF 
names(center.weights) <- centers.long$gene
center.weights.norm <- log(exp(center.weights) / sum(exp(center.weights)))
# center.weights.bin <- ifelse(center.weights.norm > -5.5, 1, 0)
center.weights.bin <- ifelse(center.weights.norm > quantile(center.weights.norm, 0.25), 1, 0)
print(head(sort(center.weights, decreasing = TRUE)))
print(head(center.weights.norm))
print(sum(center.weights.bin))
plot(density(center.weights.norm))
# plot(density(center.weights.bin))

# use weighted avg?
center.genes <- names(center.weights)
# center.genes <- sapply(names(center.weights), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)

# filter?
# jtmp <- data.frame(subset(dat.long, Gene_Name %in% top.genes & CellType == "nucleate_erythrocyte") %>% arrange(desc(zscore)))
# center.genes <- subset(centers.long, gene %in% jtmp$Gene_Name)$gene

# where are these jtmp genes 

center.dat <- data.frame(Gene_Name = center.genes, cweight = center.weights.norm, cweight.bin = center.weights.bin)

ctypes.filt <- c("T_cell", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "granulocyte", "Kit_and_Sca1-positive_hematopoietic_stem_cell")
zscores.avg.weighted <- subset(dat.long, Gene_Name %in% center.genes) %>%
  group_by(CellType) %>%
  left_join(., center.dat) %>%
  summarise(zscore = weighted.mean(zscore, w = cweight.bin)) %>%
  ungroup() %>%
  mutate(weight = exp(zscore) / sum(exp(zscore)))

ggplot(zscores.avg.weighted, aes(x = CellType, y = weight)) + geom_bar(stat = "identity") +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
  

# Redo filter top hits from only topics -----------------------------------



# try just top hits from topics
topics.ngenes <- 500
genes.filt.bytopic <- subset(terms.filt, topic %in% topics.keep & rnk <= topics.ngenes)$termgene

print(length(genes.filt.bytopic))

# further filter by marker genes?
genes.filt.bytopic.gene <- sapply(genes.filt.bytopic, function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)

filt.by.i <- which(genes.filt.bytopic.gene %in% markers)  # big majority gone, maybe zscore bad??

genes.filt.bytopic <- genes.filt.bytopic[filt.by.i]

print(length(genes.filt.bytopic))

dat.impute.sub2.bytopic <- dat.impute.sub[, genes.filt.bytopic]

# do clustering by flexmix
inmat.bytopic <- as.matrix(dat.impute.sub2.bytopic)
# log transform the inmat and then center rows

inmat.norm.bytopic <- log(inmat.bytopic * 10^6)
inmat.norm.bytopic <- sweep(inmat.norm.bytopic, MARGIN = 1, STATS = rowMeans(inmat.norm.bytopic), FUN = "-")

system.time(
  mixout.bytopic <- flexmix(inmat.bytopic ~ 1, k = nctypes, model = FLXMCmvnorm(), control = list(classify = "weighted"))
)
mixout.post.bytopic <- data.frame(cell = rownames(inmat.bytopic), as.data.frame(posterior(mixout.bytopic)))

# plot the clustering outputs
dat.clstring.bytopic <- left_join(dat.umap.mixed[[jmark]], mixout.post.bytopic)

pdf(paste0("~/data/scchic/pdfs/celltyping_analysis/celltypes_flexmix.ngenes.", topics.ngenes, ".topicsonly.", ntopicsfilt, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (i in seq(mixout.bytopic@k)){
  jcname <- paste("V", i, sep = "")
  print(jcname)
  m1 <- PlotXYWithColor(dat.clstring.bytopic, xvar = "umap1", yvar = "umap2", cname = jcname, jtitle = jcname)
  print(m1)
}
dev.off()

# this looks good??

# are NK cells marker genes?
jtop <- 22
jgenes <- subset(terms.filt, topic == jtop)$gene[1:50]

ggplot(dat.long %>%  filter(Gene_Name %in% jgenes), 
                           aes(x = CellType, y = zscore)) + geom_boxplot() + 
         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45))


# try cellassign
marker.genes.mat <- tidyr::spread(markers.keep %>% dplyr::select(Gene_Name, CellType, is.marker), key = CellType, value = is.marker) %>%
  filter(Gene_Name %in% genes.filt2) %>%
  as.data.frame()
rownames(marker.genes.mat) <- marker.genes.mat$Gene_Name
marker.genes.mat$Gene_Name <- NULL
marker.genes.mat <- as.matrix(marker.genes.mat)
sizefacs <- rep(1, nrow(inmat))

library(devtools)
dev_mode(T)
library(cellassign)



system.time(
  cas <- cellassign(exprs_obj = inmat.norm, marker_gene_info = marker.genes.mat, s = sizefacs, B = 2)
)

dev_mode(F)

# plot the hits 
cas.long <- data.frame(cell = rownames(inmat), as.data.frame(cas$mle_params$gamma), stringsAsFactors = FALSE)
cas.long <- left_join(cas.long, dat.umap.mixed[[jmark]] %>% dplyr::select(cell, umap1, umap2))

pdf("~/data/scchic/pdfs/celltyping_analysis/cellassign_markers.pdf")
for (ctype in colnames(marker.genes.mat)){
  ctype <- gsub("-", ".", ctype)
  m1 <- PlotXYWithColor(cas.long, xvar = "umap1", yvar = "umap2", cname = ctype, jtitle = ctype)
  print(m1)
}
dev.off()
