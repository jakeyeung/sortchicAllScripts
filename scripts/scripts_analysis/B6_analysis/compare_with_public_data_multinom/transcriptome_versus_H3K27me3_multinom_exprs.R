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
library(here())

library(hash)

library(scchicFuncs)

setwd(here())
source("scripts/Rfunctions/PlotFunctions.R")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- "H3K4me3"
# jmark <- "H3K4me1"
# jmark <- "H3K9me3"
jmark <- "H3K27me3"


# Constants ---------------------------------------------------------------

use.zscore <- TRUE

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
tm.result.lst <- lapply(out.objs, function(x) posterior(x$out.lda))
tm.result.lst[["H3K4me3"]] <- tm.result.stringent



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

if (use.zscore){
  dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>% as.data.frame() 
} else {
  dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-zscore), key = "celltype", value = "exprs") %>% as.data.frame() 
}
rownames(dat.norm.zscore.mat) <- dat.norm.zscore.mat$gene
dat.norm.zscore.mat$gene <- NULL


# keep only bins that are close to the TSS
terms.keep.dat <- out2.df.closest %>%
  group_by(gene) %>%
  filter(dist.to.tss == min(dist.to.tss))
terms.closest <- terms.keep.dat$region_coord

terms.sum <- terms.filt %>%
  group_by(gene) %>%
  dplyr::filter(rnk == min(rnk))

# term2gene <- hash(terms.keep.dat$region_coord, terms.keep.dat$gene)
# gene2term <- hash(terms.keep.dat$gene, terms.keep.dat$region_coord)

term2gene <- hash(terms.sum$term, terms.sum$gene)
gene2term <- hash(terms.sum$gene, terms.sum$term)

genes.keep <- rownames(dat.mat)



# Find correlated celltypes -----------------------------------------------


dat.pca <- prcomp(t(dat.mat), center = TRUE, scale. = TRUE)
dat.proj <- t(dat.mat) %*% dat.pca$rotation %*% diag(dat.pca$sdev)

plot(dat.proj[, 1], dat.proj[, 2])
text(dat.proj[, 1], dat.proj[, 2], labels = rownames(dat.proj))



# Order entropies, and check z-scores -------------------------------------

topics.mat <- tm.result.lst[[jmark]]$topics
colnames(topics.mat) <- paste("Topic", colnames(topics.mat), sep = "_")
topics.mat <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)

dat.umap.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat)

topics.sum <- OrderTopicsByEntropy(tm.result.lst[[jmark]], jquantile = 0.99) %>%
  mutate(topic = gsub("^X", "Topic_", topic))

jtopic <- "Topic_2"
ngenes.keep <- 200

terms.filt$topicstr <- paste("Topic", terms.filt$topic, sep = "_")
# paste("Topic", terms.filt$topic, sep = "_")

jtopics <- topics.sum$topic

pdf(paste0("/Users/yeung/data/scchic/pdfs/revisions/compare_", jmark, "_with_tx_quantilenorm.", Sys.Date(), ".pdf"), useDingbats = FALSE)
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

# take most interesting topics
top.nterms <- 150
topics.keep <- topics.sum$topic
# topics.keep <- c("Topic_40", "Topic_39", "Topic_27",   # monocyte
#                  "Topic_19", "Topic_8", "Topic_34", 
#                  "Topic_11", "Topic_6", "Topic_49", 
#                  "Topic_14", "Topic_1", "Topic_42")

terms.all <- unique(subset(terms.sum %>% arrange(desc(weight)), rnk <= top.nterms)$term)
genes.all <- unlist(sapply(terms.all, function(x) term2gene[[x]]))

# remove ribosomal genes
print(paste("N genes before:", length(genes.all)))
genes.all <- genes.all[!grepl("^Rp", x = genes.all)]
print(paste("N genes after:", length(genes.all)))

# genes.keep <- intersect(genes.all, rownames(dat.mat))
# terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])


# genes.all <- rownames(dat.mat)
# terms.all <- subset(terms.keep.dat, gene %in% genes.all)$region_coord
# 
# genes.keep <- sapply(terms.all, function(x) term2gene[[x]])
# terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])





# Do likelihoods ----------------------------------------------------------


# create likelihoods
jcutoff <- 1.8
plot(density(rowMeans(dat.mat)))
abline(v = jcutoff)

dat.mat.filt <- dat.mat[apply(dat.mat, 1, max) > jcutoff, ]

genes.keep <- intersect(genes.all, rownames(dat.mat.filt))
terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])

dat.mat.filt <- dat.mat.filt[genes.keep, ]

# handle zeros
zero.fill <- min(dat.mat[which(dat.mat != 0)])
dat.mat.filt[which(dat.mat.filt == 0)] <- zero.fill


# redefine genes that we keep
genes.filt <- intersect(rownames(dat.mat.filt), genes.keep)

# try using zscore
if (use.zscore){
  dat.mat.filt <- 2^(dat.norm.zscore.mat[genes.filt, ])
} else {
  dat.mat.filt <- dat.mat.filt[genes.filt, ]
}
# flip sign because represssive?


# make likelihoods
probs.lst.filt <- as.list(as.data.frame(dat.mat.filt))
# name the list just to be safe
probs.lst.filt <- lapply(probs.lst.filt, function(x){
  names(x) <- rownames(dat.mat.filt)
  return(x)
}) 



# Prepare count mat -------------------------------------------------------

count.filt <- count.dat$counts[terms.keep, ]
rownames(count.filt) <- genes.keep

# Do fits -----------------------------------------------------------------



LL.ctype.lst <- lapply(all.cells, function(cell.name){
  cell.vec <- count.filt[genes.filt, cell.name]
  # cell.vec <- cell.vec[which(cell.vec > 0)]
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    # jprob <- jprob[names(cell.vec)]
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    # return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
    # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
    # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
    return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
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
    # for repressive, take the min?
    # LL.max <- max(LL.vec)
    # p.max <- max(p.vec)
    LL.max <- min(LL.vec)
    p.max <- min(p.vec)
    best.ctype <- names(which.min(LL.vec))
  }
  dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, cell.count.downsamp = cell.count.downsamp, stringsAsFactors = FALSE)
  return(dat.tmp) 
}) %>%
  bind_rows()


pcutoff <- log(1-0.9999)

LL.dat <- LL.dat %>%
  rowwise() %>%
  mutate(ctype.stringent = ifelse(p.max <= pcutoff, ctype.pred, NA))

LL.sum <- LL.dat %>%
  group_by(ctype.stringent) %>%
  summarise(ncell = length(ctype.stringent))

print(LL.sum)

LL.dat.merge <- left_join(dat.umap.long.trajs[[jmark]], LL.dat)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = ctype.stringent)) + geom_point() +
  scale_color_manual(values = cbPalette) +
  theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.pred) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.stringent) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = LL.max/cell.size)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.pred) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Debug on a granu --------------------------------------------------------


# Test a cell -------------------------------------------------------------

plot(probs.lst.filt$granulocyte, probs.lst.filt$monocyte, pch = 20)
text(probs.lst.filt$granulocyte, probs.lst.filt$monocyte, labels = names(probs.lst.filt[[1]]))


jcell <- (subset(LL.dat.merge) %>% filter(umap2 == min(umap2)))$cell

# compare likelihoods
Lvec <- sort(LL.ctype.lst[[jcell]], decreasing = TRUE)
Pvec <- sort(p.ctype.lst[[jcell]], decreasing = TRUE)

plot(x = seq(length(Lvec)), y = Lvec)
text(x = seq(length(Lvec)), y = Lvec, labels = names(Lvec))

# plot cell on umap
ggplot(dat.umap.long.trajs[[jmark]] %>% mutate(is.cell = cell %in% jcell), aes(x = umap1, y = umap2, color = is.cell, size = is.cell)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# look at raw vector see why
cell.vec <- count.filt[genes.filt, jcell][names(probs.lst.filt[[1]])]

jgenes.remove <- c("")
par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
jctypes <- c("monocyte", "granulocyte", "megakaryocyte", "natural_killer_cell")
# jctype <- "monocyte"
for (jctype in jctypes){
  print(jctype)
  # ptmp <- probs.lst.filt[[jctype]] ^ (1 / 2)
  # ptmp <- ptmp / sum(ptmp)
  ptmp <- sqrt(probs.lst.filt[[jctype]])
  # jgenes.filt <- names(ptmp)[which(!names(ptmp) %in% jgenes.remove)]
  L <- dmultinom(x = cell.vec, prob = ptmp, log = TRUE)
  plot(x = cell.vec, y = ptmp, main = signif(L, digits = 4), pch = 20, xlab = jctype, ylab = jcell)
  text(x = cell.vec, y = ptmp, labels = names(cell.vec))
}
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)



# what is exprs of top granu genes across counts?
jsub <- count.filt[genes.filt, ]

# top 10 granu genes




# Plot top genes in UMAP --------------------------------------------------


jctype <- "megakaryocyte"

jctype <- "monocyte"
jctype <- "nucleate_erythrocyte"
jctype <- "granulocyte"

jctypes <- names(probs.lst.filt)

for (jctype in jctypes){
  jgenes.tmp <- names(sort(probs.lst.filt[[jctype]], decreasing = TRUE)[1:50])
  
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  pdf(paste0("/Users/yeung/data/scchic/pdfs/celltyping_analysis_debugging/", jmark, "_debugging_", jctype, "_jgenes_useZscore.", use.zscore, ".pdf"), useDingbats = FALSE)
  
  boxplot(dat.mat.filt[jgenes.tmp, ])
  graphics::text(seq_along(dat.mat.filt[jgenes.tmp, ]), par("usr")[3] - 0.5, labels = names(dat.mat.filt[jgenes.tmp, ]), srt = 45, adj = 1)
  
  # qplot(dat.mat.filt[jgenes.tmp, ], geom = "boxplot")
  for (jgene in jgenes.tmp){
    gene.counts <- data.frame(cell = colnames(count.filt), UMI = count.filt[jgene, ])
    LL.dat.merge.counts <- left_join(LL.dat.merge, gene.counts)
    m <- ggplot(LL.dat.merge.counts, aes(x = umap1, y = umap2, color = UMI, size = UMI)) + geom_point() +
      scale_color_viridis_c() +
      ggtitle(jgene) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  dev.off()
  
}

