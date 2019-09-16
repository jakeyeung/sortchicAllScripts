# Jake Yeung
# Date of Creation: 2019-09-15
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/transcriptome_versus_H3K4me3_multinom.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)

library(preprocessCore)
source("scripts/Rfunctions/PlotFunctions.R")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- "H3K4me3"
jmethod <- "pearson"
# Load data ---------------------------------------------------------------

# inf <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.RData"
inf <- "/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.withDist.RData"
load(inf, v=T)

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)
out.objs.stringent <- out.objs
tm.result.stringent <- posterior(out.objs$out.lda)

inf.lda.all <- "/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
load(inf.lda.all, v=T)


out.objs$H3K4me3 <- out.objs.stringent

inf.objs <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.objs, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent$H3K4me3


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  group_by(Gene_Name, CellType) %>%
  summarise(FPKM = sum(FPKM)) %>%
  rowwise() %>%
  mutate(logFPKM = log2(FPKM + 1))

# normalize across samples?
ggplot(dat.long, aes(x = CellType, y = logFPKM)) + geom_boxplot() + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
dat.mat <- tidyr::spread(dat.long %>% 
                           ungroup() %>%
                           # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>% 
                           mutate(gene = Gene_Name) %>%
                           dplyr::select(gene, CellType, logFPKM), 
                         key = CellType, value = logFPKM)  %>%
  as.data.frame()
rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL

cnames.tmp <- colnames(dat.mat)
rnames.tmp <- rownames(dat.mat)
dat.mat <- preprocessCore::normalize.quantiles(as.matrix(dat.mat), copy = TRUE)  # strong normalization,
colnames(dat.mat) <- cnames.tmp
rownames(dat.mat) <- rnames.tmp

boxplot(dat.mat)

dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

# keep only bins that are close to the TSS
terms.keep.dat <- out2.df.closest %>%
  group_by(gene) %>%
  filter(dist.to.tss == min(dist.to.tss))
terms.keep <- terms.keep.dat$region_coord

term2gene <- hash(terms.keep.dat$region_coord, terms.keep.dat$gene)
gene2term <- hash(terms.keep.dat$gene, terms.keep.dat$region_coord)

genes.keep <- rownames(dat.mat)

# Order entropies, and check z-scores -------------------------------------

topics.mat <- tm.result.stringent$topics
colnames(topics.mat) <- paste("Topic", colnames(topics.mat), sep = "_")
topics.mat <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)

dat.umap.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat)

topics.sum <- OrderTopicsByEntropy(tm.result.stringent, jquantile = 0.99) %>%
  mutate(topic = gsub("^X", "Topic_", topic))

jtopic <- "Topic_2"
ngenes.keep <- 200

terms.filt$topicstr <- paste("Topic", terms.filt$topic, sep = "_")
# paste("Topic", terms.filt$topic, sep = "_")

jtopics <- topics.sum$topic

pdf(paste0("/Users/yeung/data/scchic/pdfs/revisions/compare_", jmark, "_with_tx_quantilenorm.pdf"), useDingbats = FALSE)
for (jtopic in jtopics){
  print(jtopic)
  m.umap <- PlotXYWithColor(dat.umap.merge, xvar = "umap1", yvar = "umap2", cname = jtopic, cont.color = TRUE, jsize = 4) + scale_color_viridis_c()
  jsub <- subset(terms.filt, topicstr == jtopic & rnk <= ngenes.keep)
  print(head(jsub))
  jgenes <- jsub$gene
  m.celltype <- ggplot(subset(dat.norm.long, gene %in% jgenes),
                       aes(x = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE), y = zscore)) +
    geom_boxplot() + 
    geom_jitter(width = 0.2) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m.umap)
  print(m.celltype)
}
dev.off()




# Set up raw dataa --------------------------------------------------------


inf.raw <- paste0("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_", jmark, "_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData")
load(inf.raw, v=T)
# handle rownames
rownames(count.dat$counts) <- paste("chr", rownames(count.dat$counts), sep = "")

all.cells <- colnames(count.dat$counts)
names(all.cells) <- all.cells

genes.all <- rownames(dat.mat)
terms.all <- subset(terms.keep.dat, gene %in% genes.all)$region_coord

genes.keep <- sapply(terms.all, function(x) term2gene[[x]])
terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])

count.filt <- count.dat$counts[terms.keep, ]
rownames(count.filt) <- genes.keep

# Do likelihoods ----------------------------------------------------------


# create likelihoods
jcutoff <- 1.8
plot(density(rowMeans(dat.mat)))
abline(v = jcutoff)

dat.mat.filt <- dat.mat[rowMeans(dat.mat) > jcutoff, ]
dat.mat.filt <- dat.mat.filt[genes.keep, ]

# handle zeros
zero.fill <- min(dat.mat[which(dat.mat != 0)])
dat.mat.filt[which(dat.mat.filt == 0)] <- zero.fill


# redefine genes that we keep
genes.filt <- intersect(rownames(dat.mat.filt), genes.keep)

dat.mat.filt <- dat.mat.filt[genes.filt, ]

# make likelihoods
probs.lst.filt <- as.list(as.data.frame(dat.mat.filt))
# name the list just to be safe
probs.lst.filt <- lapply(probs.lst.filt, function(x){
  names(x) <- rownames(dat.mat.filt)
  return(x)
}) 



# Do fits -----------------------------------------------------------------



LL.ctype.lst <- lapply(all.cells, function(cell.name){
  cell.vec <- count.filt[genes.filt, cell.name]
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
  })
})




# Summarize fits ----------------------------------------------------------

# calculate probability of model given data 
p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)

cell.counts <- Matrix::colSums(count.dat$counts) / 5
cell.counts.downsamp <- Matrix::colSums(count.dat$counts) / 5

cell.names <- names(all.cells)

# summaize
LL.dat <- lapply(cell.names, function(cname){
  LL.vec <- LL.ctype.lst[[cname]]
  p.vec <- p.ctype.lst[[cname]]
  cell.count = cell.counts[[cname]]
  cell.count.downsamp = cell.counts.downsamp[[cname]]
  if (all(is.infinite(LL.vec))){
    LL.max <- NA
    p.max <- NA
    best.ctype <- NA
  } else {
    LL.max <- max(LL.vec)
    p.max <- max(p.vec)
    best.ctype <- names(which.max(LL.vec))
  }
  dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, cell.count.downsamp = cell.count.downsamp, stringsAsFactors = FALSE)
  return(dat.tmp) 
}) %>%
  bind_rows()

LL.sum <- LL.dat %>%
  group_by(ctype.pred) %>%
  summarise(ncell = length(ctype.pred))

print(LL.sum)

LL.dat.merge <- left_join(dat.umap.long.trajs.stringent[[jmark]], LL.dat)


ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.pred) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# why monocytes?
granu.genes <- subset(terms.filt, topic == "42")$gene[1:10]

probs.lst.filt$granulocyte[granu.genes]
probs.lst.filt$monocyte[granu.genes]

# why does top 10 granu genes show separation ? 
ggplot(subset(dat.norm.long, gene %in% granu.genes), 
       aes(x = forcats::fct_reorder(celltype, zscore, median), y = zscore, group = celltype)) + 
  geom_point() +  geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

