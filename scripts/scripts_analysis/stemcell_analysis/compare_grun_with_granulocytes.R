# Jake Yeung
# Date of Creation: 2019-11-24
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/compare_grun_with_granulocytes.R
# do it all 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(hash)
library(igraph)
library(umap)
library(Seurat)
library(scchicFuncs)

library(preprocessCore)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(ggrepel)

library(DESeq2)

library(Matrix)

jmark <- "H3K4me3"

LoadTSSCounts <- function(inf.tss){
  mat <- fread(inf.tss)
  colnames(mat)[1:4] <- unlist(mat[1, 1:4])
  mat <- mat[-1, ]
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$bname
  return(mat)
}


# Load trajs --------------------------------------------------------------

inf.objs <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.objs, v=T)



# Load grun  --------------------------------------------------------------


# load dat ----------------------------------------------------------------

inf.clst <- "/Users/yeung/data/scchic/public_data/bone_marrow_grun/clustering_bone_marrow.txt.gz"
meta <- fread(inf.clst, header = FALSE, col.names = c("cell", "cluster")) %>%
  as.data.frame()
rownames(meta) <- meta$cell

# inf <- "/Users/yeung/data/scchic/public_data/bone_marrow_grun/transcript_counts_normalized_bone_marrow.txt.gz"
inf <- "/Users/yeung/data/scchic/public_data/bone_marrow_grun/GSE76983_expdata_BMJhscC.csv.gz"
dat <- fread(inf, header = TRUE) %>%
  as.data.frame()
gene.names <- dat$GENEID
# gene.names <- sapply(gene.names, function(x) strsplit(x, "__")[[1]][[1]])

# inf.test <- "/Users/yeung/data/scchic/public_data/bone_marrow_grun/transcript_counts_bone_marrow.txt.gz"
inf.test <- "/Users/yeung/data/scchic/public_data/bone_marrow_grun/transcript_counts_bone_marrow.xls"
dat.test <- fread(inf.test, header = FALSE) %>%
  as.data.frame()
rownames(dat.test) <- rownames(dat)

cols.keep <- which(colnames(dat) %in% meta$cell)
dat.filt <- dat[, c(1, cols.keep)]



rownames(dat.filt) <- gene.names
dat.filt$GENEID <- NULL

# guess
colnames(dat.test) <- colnames(dat.filt)
rownames(dat.test) <- rownames(dat.filt)

identical(dat.filt, dat.test)

# jdiff <- unlist(dat.filt - dat.test)





# Create averages ---------------------------------------------------------

count.long <- data.table::melt(data.frame(gene = rownames(dat.filt), dat.filt, stringsAsFactors = FALSE), variable.name = "cell", value.name = "counts")


clstrs.keep <- c(1, 11, 3, 2, 14, 12)

# rename clusters
clstrs.rename <- c("stem", "neutro1", "neutro2", "neutro3", "neutro4", "neutro5")

clstrs.hash <- hash(clstrs.keep, clstrs.rename)
clstrs.hash[["Other"]] <- "Other"

# average across  clusters
meta.neutro <- meta %>% mutate(is.neutro = cluster %in% clstrs.keep, neutro.stage = ifelse(is.neutro, cluster, "Other"))
meta.neutro$clstr.name <- sapply(as.character(meta.neutro$neutro.stage), function(x) clstrs.hash[[x]])
rownames(meta.neutro) <- meta.neutro$cell

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# sum across cells
count.long.merge <- left_join(count.long, subset(meta.neutro, select = c(-is.neutro, -neutro.stage))) %>%
  group_by(clstr.name, gene) %>%
  summarise(count.sum = sum(counts)) %>%
  group_by(clstr.name) %>%
  mutate(count.norm = 10^6 * count.sum / sum(count.sum)) %>%
  rowwise() %>%
  mutate(logcount.norm = log2(count.norm + 1)) %>%
  group_by(gene) %>%
  mutate(zscore = scale(logcount.norm, center = TRUE, scale = TRUE))

count.mat.bulk <- count.long.merge %>%
  dplyr::select(gene, clstr.name, count.sum) %>%
  tidyr::spread(key = clstr.name, value = count.sum) %>%
  as.data.frame()


rownames(count.mat.bulk) <- count.mat.bulk$gene
count.mat.bulk$gene <- NULL

boxplot(count.mat.bulk, main = "Before quant norm")

# remove outliers, low number of genes?
cname.remove <- "neutro4"
# cname.remove <- ""
cname.keep.i <- which(!colnames(count.mat.bulk) %in% cname.remove)
count.mat.bulk <- count.mat.bulk[, cname.keep.i]


# do DESeq2 variance stabilization 

metadata <- data.frame(ctype = colnames(count.mat.bulk), stringsAsFactors = FALSE)
rownames(metadata) <- colnames(count.mat.bulk)

count.mat.bulk.int <- matrix(as.integer(as.matrix(count.mat.bulk)), 
                             nrow = nrow(count.mat.bulk), ncol = ncol(count.mat.bulk), 
                             dimnames = list(rownames(count.mat.bulk), colnames(count.mat.bulk)))
dds <- DESeqDataSetFromMatrix(countData = count.mat.bulk.int, colData = metadata, design = ~1)
vsd <- vst(dds)

boxplot(assay(vsd))

exprs.cutoff <- 3
plot(density(assay(vsd)))
abline(v = exprs.cutoff)

# remove low counts
dat.vst.filt <- assay(vsd)[which(rowMeans(assay(vsd)) > exprs.cutoff), ]

rnames.orig <- rownames(dat.vst.filt)
cnames.orig <- colnames(dat.vst.filt)

# count.mat.bulk <- preprocessCore::normalize.quantiles(as.matrix(assay(vsd)), copy = TRUE)
count.mat.bulk.norm <- preprocessCore::normalize.quantiles(as.matrix(dat.vst.filt), copy = TRUE)

rownames(count.mat.bulk.norm) <- rnames.orig
colnames(count.mat.bulk.norm) <- cnames.orig

# split chromosome location from gene
rownames(count.mat.bulk.norm) <- sapply(rownames(count.mat.bulk.norm), function(x) strsplit(x, "__")[[1]][[1]])

boxplot(count.mat.bulk.norm, main = "After quant norm")

dat.pca <- prcomp(t(count.mat.bulk.norm), center = TRUE, scale. = TRUE)

dat.pca.proj <- as.data.frame(t(count.mat.bulk.norm) %*% dat.pca$rotation %*% diag(dat.pca$sdev))
colnames(dat.pca.proj) <- paste("PC", seq(ncol(dat.pca.proj)), sep = "")
dat.pca.proj$clstr.name <- rownames(dat.pca.proj)

ggplot(dat.pca.proj, aes(x = PC1, y = PC2, label = clstr.name)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_text_repel()



# Set up counts -----------------------------------------------------------

# # remove bad genes
bad.genes <- "Gm1821"
# count.mat.bulk.norm <- count.mat.bulk.norm[!grepl(bad.genes, rownames(count.mat.bulk.norm)), ]
# rownames(count.mat.bulk.norm) <- sapply(rownames(count.mat.bulk.norm), function(x) strsplit(x, "--")[[1]][[1]])

# make long
count.mat.bulk.norm.long <- melt(data.frame(gene = rownames(count.mat.bulk.norm), count.mat.bulk.norm, stringsAsFactors = FALSE), variable.name = "clstr.name", value.name = "counts") %>%
  group_by(gene) %>%
  mutate(zscore = scale(counts, center = TRUE, scale = TRUE))

exprs.mat <- count.mat.bulk.norm.long %>%
  dplyr::select(gene, clstr.name, counts) %>%
  tidyr::spread(key = clstr.name, value = counts) %>%
  as.data.frame()
rownames(exprs.mat) <- exprs.mat$gene
exprs.mat <- exprs.mat[!grepl(bad.genes, rownames(exprs.mat)), ]
exprs.mat$gene <- NULL

zscore.mat <- count.mat.bulk.norm.long %>%
  dplyr::select(gene, clstr.name, zscore) %>%
  tidyr::spread(key = clstr.name, value = zscore) %>%
  as.data.frame()
rownames(zscore.mat) <- zscore.mat$gene
zscore.mat <- zscore.mat[!grepl(bad.genes, rownames(zscore.mat)), ]
zscore.mat$gene <- NULL
# rownames(zscore.mat) <- sapply(rownames(zscore.mat), function(x) strsplit(x, "--")[[1]][[1]])

# dat.mat <- zscore.mat


# Load raw data  ----------------------------------------------------------

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)
out.objs.stringent <- out.objs
# tm.result.stringent <- posterior(out.objs$out.lda)

infs.tss <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/countTables_TSS/countTables_TSS", pattern = "*H3K4me3.*.csv", full.names = TRUE)



mats.tss <- lapply(infs.tss, LoadTSSCounts)

rows.common <- Reduce(intersect, lapply(mats.tss, rownames))
cells.keep <- colnames(out.objs.stringent$count.mat)

mats.tss.filt <- lapply(mats.tss, function(x){
  rows.i <- which(rownames(x) %in% rows.common)
  cols.i <- which(colnames(x) %in% cells.keep)
  return(x[rows.i, cols.i])
}) %>%
  do.call(cbind, .)

mats.tss.filt <- Matrix(as.matrix(mats.tss.filt), sparse = TRUE)
mats.tss.filt[which(is.na(mats.tss.filt))] <- 0

# genes.chic <- sapply(rownames(mats.tss.filt), function(x) strsplit(x, ";")[[1]][[2]])
mats.tss.filt.collapsed <- CollapseRowsByGene(mats.tss.filt, as.long = FALSE)



# Do likelihoods ----------------------------------------------------------

genes.all <- rownames(mats.tss.filt.collapsed)
genes.keep <- intersect(genes.all, rownames(zscore.mat))

# dat.mat.filt <- zscore.mat[genes.keep, ]
raw.mat.filt <- mats.tss.filt.collapsed[genes.keep, ]

# exponentiate
# dat.mat.filt <- 2^zscore.mat[genes.keep, ]

dat.mat.filt <- exprs.mat[genes.keep, ]

# try using zscore
# dat.mat.filt <- 2^zscore.mat[genes.keep, cname.keep.i]  # exp is worse than 2

# dat.mat.filt <- dat.norm.zscore.mat[genes.filt, ]
# renormalize so no zeros??
# dat.mat.filt <- sweep(dat.mat.filt, MARGIN = 2, STATS = apply(dat.mat.filt, 2, min), FUN = "-")
# dat.mat.filt <- as.matrix(dat.mat.filt)

# handle zeros
zero.fill <- min(as.matrix(dat.mat.filt)[which(as.matrix(dat.mat.filt) != 0)])
dat.mat.filt[which(dat.mat.filt == 0)] <- zero.fill

# make likelihoods
probs.lst.filt <- as.list(as.data.frame(dat.mat.filt))
# name the list just to be safe
probs.lst.filt <- lapply(probs.lst.filt, function(x){
  names(x) <- rownames(dat.mat.filt)
  return(x)
})

plot(probs.lst.filt$neutro5, probs.lst.filt$neutro3, pch = 20)
text(probs.lst.filt$neutro5, probs.lst.filt$neutro3, labels = names(probs.lst.filt$neutro3))


# Fit raw data ------------------------------------------------------------

all.cells <- colnames(raw.mat.filt)
names(all.cells) <- all.cells

jexppower <- 0.5
LL.ctype.lst <- FitMultinoms(raw.mat.filt, all.cells, probs.lst.filt, exppower = jexppower)
LL.dat <- SummarizeMultinomFits(LL.ctype.lst, raw.mat.filt, all.cells) %>%
  mutate(is.stem = FALSE)

dat.umap.long.merge <- left_join(dat.umap.long, LL.dat, by = "cell")


ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = ctype.pred)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


dat.umap.long.merge.traj <- left_join(dat.umap.long.merge, dat.umap.long.trajs[[jmark]])

ggplot(dat.umap.long.merge.traj %>% filter(granu & p.max > log(0.9)), aes(x = umap1, y = umap2, color = ctype.pred)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + facet_wrap(~ctype.pred)

ggplot(dat.umap.long.merge.traj %>% filter(granu & p.max > log(0.9)), aes(x = umap1, y = umap2, color = ctype.pred)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

# plot some markers from Grun


jgenes <- c("Ngp", "Retnlg", "Elane", "Junb", "S100a9", "Top2a", "Tuba1b", "Rpl7a", "Cdk1")

marker.exprs <- as.data.frame(t(raw.mat.filt[jgenes, ]))

dat.umap.long.merge.traj.marker <- left_join(dat.umap.long.merge.traj, data.frame(cell = rownames(marker.exprs), marker.exprs, stringsAsFactors = FALSE))

PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Cdk1", jsize = 3) + scale_color_viridis_c()

PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Rpl7a", jsize = 3) + scale_color_viridis_c()

PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "S100a9", jsize = 3) + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Junb", jsize = 3) + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Top2a", jsize = 3) + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Tuba1b", jsize = 3) + scale_color_viridis_c()


PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Retnlg") + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Elane") + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Ngp") + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% filter(granu), xvar = "umap1", yvar = "umap2", cname = "Ngp") + scale_color_viridis_c()

PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette)

# 5 versus 1 versus 8?
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% rowwise() %>% mutate(louvain = ifelse(louvain == 5, TRUE, FALSE)), xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette)
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% rowwise() %>% mutate(louvain = ifelse(louvain == 1, TRUE, FALSE)), xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette)
PlotXYWithColor(dat.umap.long.merge.traj.marker %>% rowwise() %>% mutate(louvain = ifelse(louvain == 8, TRUE, FALSE)), xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette)


# Get top 5 versus all ----------------------------------------------------

x <- probs.lst.filt$neutro5
y <- rowMeans(data.frame(probs.lst.filt$neutro1, probs.lst.filt$neutro2, probs.lst.filt$neutro3))
jgenes <- names(which((x + y) / 2 > 6 & (x - y) > 2))

x <- probs.lst.filt$neutro3
y <- rowMeans(data.frame(probs.lst.filt$neutro1, probs.lst.filt$neutro2, probs.lst.filt$neutro5))
jgenes <- names(which((x + y) / 2 > 4 & (x - y) > 2))

x <- probs.lst.filt$neutro2
y <- rowMeans(data.frame(probs.lst.filt$neutro1, probs.lst.filt$neutro3, probs.lst.filt$neutro5))
jgenes <- names(which((x + y) / 2 > 4 & (x - y) > 1.5))

x <- probs.lst.filt$neutro1
y <- rowMeans(data.frame(probs.lst.filt$neutro1, probs.lst.filt$neutro2, probs.lst.filt$neutro5))
jgenes <- names(which((x + y) / 2 > 6 & (x - y) > 1))

# x <- probs.lst.filt$stem
# y <- rowMeans(data.frame(probs.lst.filt$neutro1, probs.lst.filt$neutro2, probs.lst.filt$neutro3, probs.lst.filt$neutro5))

x <- probs.lst.filt$Other
y <- rowMeans(data.frame(probs.lst.filt$neutro1, probs.lst.filt$neutro2, probs.lst.filt$neutro5))
jgenes <- names(which((x + y) / 2 > 8 & (x - y) < -2))

jgenes <- c("S100a8", "S100a9", "S100a11")

# plot(x, y, pch = 20)
# text(x, y, labels = names(x))

# MA
plot((x + y) / 2, (x - y), pch = 20)
xfilt <- x[names(x) %in% jgenes]
yfilt <- y[names(y) %in% jgenes]
text((xfilt + yfilt) / 2, (xfilt - yfilt), labels = names(xfilt))
abline(h = 0, col = "blue")

# take top hits
# jgenes <- names(which((x + y) / 2 > 8 & (x - y) > 3))

marker.exprs <- data.frame(cell = colnames(raw.mat.filt), clstr5exprs = rowSums(as.data.frame(t(raw.mat.filt[jgenes, ]))) / colSums(raw.mat.filt), stringsAsFactors = FALSE)
marker.exprs <- data.frame(cell = colnames(raw.mat.filt), clstr5exprs = log10(10^6 * rowSums(as.data.frame(t(raw.mat.filt[jgenes, ]))) / colSums(raw.mat.filt)) + 1, stringsAsFactors = FALSE)
dat.umap.long.merge.traj.marker <- left_join(dat.umap.long.merge.traj, data.frame(cell = rownames(marker.exprs), marker.exprs, stringsAsFactors = FALSE))

PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstr5exprs", jsize = 3) + 
  scale_color_viridis_c()


# Plot for each cluster  ----------------------------------------------------------

neutro.clsts <- grep("^neutro", names(probs.lst.filt), value = TRUE)

PlotClstrGenes <- function(jclstname, min.mean = 6, min.diff = 2, plot.MA.only = FALSE, less.than = FALSE){
  jclstother.i <- which(!neutro.clsts %in% jclstname)
  x <- probs.lst.filt[[jclstname]]
  y <- rowMeans(data.frame(probs.lst.filt[jclstother.i]))
  
  plot((x + y) / 2, (x - y), pch = 20, xlab = c("Mean exprs across neutros in pseudobulk scRNAseq"), ylab = paste("Log fold change between", jclstname, "vs others"))
  if (plot.MA.only){
    return(NULL)
  }
  if (!less.than){
    jgenes <- names(which((x + y) / 2 > min.mean & (x - y) > min.diff))  # manual 
  } else {
    jgenes <- names(which((x + y) / 2 > min.mean & (x - y) < min.diff))  # manual 
  }
  xfilt <- x[names(x) %in% jgenes]
  yfilt <- y[names(y) %in% jgenes]
  text((xfilt + yfilt) / 2, (xfilt - yfilt), labels = names(xfilt))
  abline(h = 0, col = "blue")
  
  marker.exprs <- data.frame(cell = colnames(raw.mat.filt), clstrexprs = rowSums(as.data.frame(t(raw.mat.filt[jgenes, ]))) / colSums(raw.mat.filt), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(clstrexprs.log = log10(clstrexprs + 1))
  dat.umap.long.merge.traj.marker <- left_join(dat.umap.long.merge.traj, marker.exprs, stringsAsFactors = FALSE)
  m <- PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs", jsize = 3, jtitle = jclstname) + scale_color_viridis_c()
  print(m)
  m <- PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs.log", jsize = 3, jtitle = jclstname) + scale_color_viridis_c()
  print(m)
  return(dat.umap.long.merge.traj.marker)
}

# clst5
pdf("/Users/yeung/data/scchic/pdfs/stemcell_analysis.redo/integrate_with_grun.pdf", useDingbats = FALSE)
jclstname <- "neutro5"
dat.umap.long.merge.traj.marker <- PlotClstrGenes(jclstname, min.mean = 0, min.diff = 1, plot.MA.only = FALSE)

jclstname <- "neutro3"
dat.umap.long.merge.traj.marker <- PlotClstrGenes(jclstname, min.mean = 0, min.diff = 1, plot.MA.only = FALSE)

jclstname <- "neutro2"
dat.umap.long.merge.traj.marker <- PlotClstrGenes(jclstname, min.mean = 0, min.diff = 1, plot.MA.only = FALSE)

jclstname <- "neutro1"
dat.umap.long.merge.traj.marker <- PlotClstrGenes(jclstname, min.mean = 0, min.diff = 1, plot.MA.only = FALSE)
dev.off()

# Find middle genes -------------------------------------------------------

jlouv <- 1  # middle cluster
jlouv.other <- c(5, 8)  # termianting and beginning cluster, respectively

# jlouv <- 5  # middle cluster
# jlouv.other <- c(1, 8)  # termianting and beginning cluster, respectively 
cells.filt <- subset(dat.umap.long.merge.traj.marker, louvain == jlouv)$cell
cells.filt.other <- subset(dat.umap.long.merge.traj.marker, louvain %in% jlouv.other)$cell

x <- log2(1000 * rowSums(raw.mat.filt[, cells.filt]) / sum(raw.mat.filt[, cells.filt]) + 1)
y <- log2(1000 * rowSums(raw.mat.filt[, cells.filt.other])  / sum(raw.mat.filt[, cells.filt.other]) + 1)

# plot(mat.filt, mat.filt.other, pch = 20) 
# text(mat.filt, mat.filt.other, labels = names(mat.filt), pch = 20) 

plot( (x + y) / 2, (x - y) , pch = 20)
# jgenes <- names(which((x + y) / 2 > 0.35 & abs((x - y)) > 0.1))  # manual 
jgenes <- names(which((x + y) / 2 > 0.3 & (x - y) > 0.1))  # manual 
xfilt <- x[names(x) %in% jgenes]
yfilt <- y[names(y) %in% jgenes]
text((xfilt + yfilt) / 2, (xfilt - yfilt), labels = names(xfilt))

marker.exprs <- data.frame(cell = colnames(raw.mat.filt), clstrexprs = rowSums(as.data.frame(t(raw.mat.filt[jgenes, ]))) / colSums(raw.mat.filt), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(clstrexprs.log = log10(clstrexprs + 1))
dat.umap.long.merge.traj.marker <- left_join(dat.umap.long.merge.traj, marker.exprs, stringsAsFactors = FALSE)
m <- PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs", jsize = 3, jtitle = c("Genes high in middle")) + scale_color_viridis_c()
print(m)

# where are these genes in the grun dataset?
boxplot(exprs.mat[jgenes, ])


# 
# 
# jclstother.i <- which(!neutro.clsts %in% jclstname)
# x <- probs.lst.filt[[jclstname]]
# y <- rowMeans(data.frame(probs.lst.filt[jclstother.i]))
# jgenes <- names(which((x + y) / 2 > 6 & (x - y) > 2))  # manual 
# 
# plot((x + y) / 2, (x - y), pch = 20, xlab = c("Mean exprs across neutros in pseudobulk scRNAseq"), ylab = paste("Log fold change between", jclstname, "vs others"))
# xfilt <- x[names(x) %in% jgenes]
# yfilt <- y[names(y) %in% jgenes]
# text((xfilt + yfilt) / 2, (xfilt - yfilt), labels = names(xfilt))
# abline(h = 0, col = "blue")
# 
# marker.exprs <- data.frame(cell = colnames(raw.mat.filt), clstrexprs = rowSums(as.data.frame(t(raw.mat.filt[jgenes, ]))) / colSums(raw.mat.filt), stringsAsFactors = FALSE) %>%
#   rowwise() %>%
#   mutate(clstrexprs.log = log10(clstrexprs + 1))
# dat.umap.long.merge.traj.marker <- left_join(dat.umap.long.merge.traj, marker.exprs, stringsAsFactors = FALSE)
# PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs", jsize = 3) + scale_color_viridis_c()
# PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs.log", jsize = 3) + scale_color_viridis_c()
# 
# 
# # clst3
# jclstname <- "neutro3"
# jclstother.i <- which(!neutro.clsts %in% jclstname)
# x <- probs.lst.filt[[jclstname]]
# y <- rowMeans(data.frame(probs.lst.filt[jclstother.i]))
# jgenes <- names(which((x + y) / 2 > 4 & (x - y) > 2))
# 
# plot((x + y) / 2, (x - y), pch = 20, xlab = c("Mean exprs across neutros in pseudobulk scRNAseq"), ylab = paste("Log fold change between", jclstname, "vs others"))
# xfilt <- x[names(x) %in% jgenes]
# yfilt <- y[names(y) %in% jgenes]
# text((xfilt + yfilt) / 2, (xfilt - yfilt), labels = names(xfilt))
# abline(h = 0, col = "blue")
# 
# marker.exprs <- data.frame(cell = colnames(raw.mat.filt), clstrexprs = rowSums(as.data.frame(t(raw.mat.filt[jgenes, ]))) / colSums(raw.mat.filt), stringsAsFactors = FALSE) %>%
#   rowwise() %>%
#   mutate(clstrexprs.log = log10(clstrexprs + 1))
# dat.umap.long.merge.traj.marker <- left_join(dat.umap.long.merge.traj, marker.exprs, stringsAsFactors = FALSE)
# PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs", jsize = 3) + scale_color_viridis_c()
# PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs.log", jsize = 3) + scale_color_viridis_c()
# 
# 
# 
# 
# # PlotXYWithColor(dat.umap.long.merge.traj.marker, xvar = "umap1", yvar = "umap2", cname = "clstrexprs.log", jsize = 3) + scale_color_viridis_c()
# 
# 
# 
# # Do more -----------------------------------------------------------------
# 
# 
# 
# # take top genes??
# 



# correlate the raw counts?
jlouv <- 5
cells.filt <- subset(dat.umap.long.merge.traj.marker, louvain == jlouv)$cell

par(mfrow=c(2, 4), pty = "s")
for (jclstr in c("neutro1", "neutro2", "neutro3", "neutro4", "neutro5", "stem", "Other")){
  if (jclstr %in% cname.remove){
    print("Skipping")
    print(jclstr)
    next
  }
  # plot(rowSums(raw.mat.filt[genes.keep, cells.filt]), zscore.mat[genes.keep, jclstr], pch = 20, main = jclstr)
  print(jclstr)
  x <- exprs.mat[genes.keep, jclstr]
  y <- rowSums(raw.mat.filt[genes.keep, cells.filt])
  jcor <- cor(x, y, method = "spearman")
  plot(x, y, pch = 20, main = paste(jclstr, "vs", jlouv, "\n", signif(jcor, digits = 2)))
}
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)



# Find differential genes across each cluster -----------------------------



# Check whether these genes associate with any neutro clusters  -----------






# correlate?


# plot(rowSums(raw.mat.filt[genes.keep, cells.filt]), zscore.mat[genes.keep, "neutro5"], pch = 20)
# 
# 
# # Summarize fits ----------------------------------------------------------
# 
# 
# # calculate probability of model given data 
# p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)
# 
# cell.counts <- Matrix::colSums(out.objs$count.mat) / 5
# cell.counts.downsamp <- Matrix::colSums(out.objs$count.mat) / 5
# 
# cell.names <- names(all.cells)
# 
# # summaize
# LL.dat <- lapply(cell.names, function(cname){
#   LL.vec <- LL.ctype.lst[[cname]]
#   p.vec <- p.ctype.lst[[cname]]
#   cell.count = cell.counts[[cname]]
#   cell.count.downsamp = cell.counts.downsamp[[cname]]
#   if (all(is.infinite(LL.vec))){
#     LL.max <- NA
#     p.max <- NA
#     best.ctype <- NA
#   } else {
#     LL.max <- max(LL.vec)
#     p.max <- max(p.vec)
#     best.ctype <- names(which.max(LL.vec))
#   }
#   dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, cell.count.downsamp = cell.count.downsamp, stringsAsFactors = FALSE)
#   return(dat.tmp) 
# }) %>%
#   bind_rows()
# 
# # be stringent with the predictions 
# 
# 
# LL.sum <- LL.dat %>%
#   group_by(ctype.pred) %>%
#   summarise(ncell = length(ctype.pred))
# 
# print(LL.sum)
# 
# p.filt <- log(0.99)
# LL.dat <- LL.dat %>%
#   rowwise() %>%
#   mutate(ctype.stringent = ifelse(p.max >= p.filt, ctype.pred, NA))
# 
# LL.dat.merge <- left_join(dat.umap.pred.merged , LL.dat)
# 
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# m.umap.celltype <- ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = ctype.stringent)) + geom_point() + 
#   scale_color_manual(values = cbPalette, na.value = "grey95") +
#   theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# print(m.umap.celltype)
# 
# 
