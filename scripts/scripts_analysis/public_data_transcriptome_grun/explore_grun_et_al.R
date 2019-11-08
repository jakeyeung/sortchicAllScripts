# Jake Yeung
# Date of Creation: 2019-10-02
# File: ~/projects/scchic/scripts/scripts_analysis/public_data_transcriptome_grun/explore_grun_et_al.R
# Explore grun

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

# Functions ---------------------------------------------------------------



# Functions ---------------------------------------------------------------

DoFishersTestCtype <- function(jctype, celltype.fc){
  conting.tab <- spread(celltype.fc %>%
                          mutate(ctype.pred = ifelse(ctype.pred == jctype, jctype, "zNot")) %>%
                          dplyr::select(is.stem, ctype.pred, cell.count) %>%
                          group_by(ctype.pred, is.stem) %>%
                          summarise(cell.count = sum(cell.count)),
                        key = is.stem, value = cell.count) %>%
    as.data.frame() %>%
    ungroup()
  if (any(is.na(conting.tab))){
    warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
    return(NULL)
  }
  # rownames
  rownames(conting.tab) <- conting.tab$ctype.pred
  conting.tab$ctype.pred <- NULL
  
  # make odds ratio interpretable:( zNot_treat / ctype_treat ) / ( zNot_ctrl / ctype_ctrl )
  conting.tab <- t(conting.tab)
  # swap rows
  conting.tab <- conting.tab[c(2, 1), ]
  print(conting.tab)
  hyp.test <- fisher.test(x = conting.tab)
  return(hyp.test)
}


DoFishersTestCtypeStringent <- function(jctype, celltype.fc){
  conting.tab <- spread(celltype.fc %>%
                          mutate(ctype.stringent = ifelse(ctype.stringent == jctype, jctype, "zNot")) %>%
                          dplyr::select(is.stem, ctype.stringent, cell.count) %>%
                          group_by(ctype.stringent, is.stem) %>%
                          summarise(cell.count = sum(cell.count)),
                        key = is.stem, value = cell.count) %>%
    as.data.frame() %>%
    ungroup()
  if (any(is.na(conting.tab))){
    warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
    return(NULL)
  }
  # rownames
  rownames(conting.tab) <- conting.tab$ctype.stringent
  conting.tab$ctype.stringent <- NULL
  
  # make odds ratio interpretable:( zNot_treat / ctype_treat ) / ( zNot_ctrl / ctype_ctrl )
  conting.tab <- t(conting.tab)
  # swap rows
  conting.tab <- conting.tab[c(2, 1), ]
  print(conting.tab)
  hyp.test <- fisher.test(x = conting.tab)
  return(hyp.test)
}


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



# Do analysis seurat ------------------------------------------------------

bm <- CreateSeuratObject(dat.filt, project = "BM", assay = "RNA", meta.data = meta)

bm <- SCTransform(bm)

bm <- RunPCA(bm, verbose = FALSE)
bm <- RunUMAP(bm, dims = 1:30, verbose = FALSE)

bm <- FindNeighbors(bm, dims = 1:30, verbose = FALSE)
bm <- FindClusters(bm, verbose = FALSE)

DimPlot(bm, label = TRUE, group.by = "cluster") + scale_color_viridis_d()

FeaturePlot(bm, features = c("Elane--chr10", "Retnlg--chr16", "Ngp--chr9"), pt.size = 0.2, ncol = 3)


meta <- left_join(meta, data.frame(cell = rownames(bm@reductions$umap@cell.embeddings), bm@reductions$umap@cell.embeddings, stringsAsFactors = FALSE))

# Downstream --------------------------------------------------------------

# plot over pseudotime? 
count.long.log <- as.data.frame(bm@assays$SCT@data) %>%
  ungroup() %>%
  mutate(gene = rownames(bm@assays$SCT@data)) %>%
  tidyr::gather(key = "cell", value = "logcount.norm", -gene)

count.long.lin <- as.data.frame(bm@assays$SCT@counts) %>%
  ungroup() %>%
  mutate(gene = rownames(bm@assays$SCT@counts)) %>%
  tidyr::gather(key = "cell", value = "count.norm", -gene)

count.long.lin.raw <- as.data.frame(bm@assays$RNA@counts) %>%
  ungroup() %>%
  mutate(gene = rownames(as.data.frame(bm@assays$RNA@counts))) %>%
  tidyr::gather(key = "cell", value = "count.raw", -gene)

count.long <- left_join(left_join(count.long.log, count.long.lin), count.long.lin.raw)


# keep neutrophils only
clstrs.keep <- c(1, 11, 3, 2, 14, 12)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(meta %>% mutate(is.neutro = cluster %in% clstrs.keep, neutro.stage = ifelse(is.neutro, cluster, NA)), 
       aes(x = UMAP_1, y = UMAP_2, color = as.character(neutro.stage))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette, na.value = "grey85")


# Define genes ------------------------------------------------------------


# define genes.all

inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedAll_2019-09-29.CountThres0.K-30_35_50.OutObjs.RData"
load(inf.lda, v=T)
annot.out <- AnnotateBins(terms.mat = out.objs$tm.result$terms)

terms.sum <- annot.out$terms.filt %>%
  group_by(gene) %>%
  dplyr::filter(rnk == min(rnk))

term2gene <- hash(terms.sum$term, terms.sum$gene)
gene2term <- hash(terms.sum$gene, terms.sum$term)

genes.keep <- rownames(dat.mat)

top.nterms <- 150
terms.all <- unique(subset(terms.sum %>% arrange(desc(weight)), rnk <= top.nterms)$term)
genes.all <- unlist(sapply(terms.all, function(x) term2gene[[x]]))


# Set up tables -----------------------------------------------------------


# plot out
jsettings <- umap.defaults
jsettings$n_neighbors <- 25
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
# umap.out <- umap(lsi.out$u, config = jsettings)

inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.Robj"
load(inf.lda, v=T)

out.lda <- out.lda[[3]]

umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
dat.umap.long.lda <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

dat.umap.long.lda$is.stem <- sapply(dat.umap.long.lda$cell, function(x) grepl("stem-cell", x))

ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2)) + geom_point() +
  facet_wrap(~is.stem) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.long.lda <- DoLouvain(topics.mat = posterior(out.lda)$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.lda)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.louvain <- ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() +
  facet_wrap(~is.stem) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


inf.proj <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/projections/H3K4me1_matsMergedNonenriched.matsMergedEnriched.RData"
load(inf.proj, v=T)

dat.pred <- predict(umap.out, out.lda.predict$topics)

dat.pred.long <- data.frame(cell = rownames(dat.pred), umap1 = dat.pred[, 1], umap2 = dat.pred[, 2], stringsAsFactors = FALSE) %>%
  mutate(is.stem = TRUE)

dat.umap.pred.merged <- bind_rows(subset(dat.umap.long.lda, select = -louvain), dat.pred.long)

ggplot(dat.umap.pred.merged, aes(x = umap1, y = umap2)) + geom_point() +
  facet_wrap(~is.stem) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Create averages ---------------------------------------------------------


# rename clusters
clstrs.rename <- c("stem", "neutro1", "neutro2", "neutro3", "neutro4", "neutro5")

clstrs.hash <- hash(clstrs.keep, clstrs.rename)
clstrs.hash[["Other"]] <- "Other"

# average across  clusters
meta.neutro <- meta %>% mutate(is.neutro = cluster %in% clstrs.keep, neutro.stage = ifelse(is.neutro, cluster, "Other"))
meta.neutro$clstr.name <- sapply(as.character(meta.neutro$neutro.stage), function(x) clstrs.hash[[x]])
rownames(meta.neutro) <- meta.neutro$cell

ggplot(meta.neutro, 
       aes(x = UMAP_1, y = UMAP_2, color = clstr.name)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette, na.value = "grey85")

# sum across cells
count.long.merge <- left_join(count.long, subset(meta.neutro, select = c(-is.neutro, -neutro.stage))) %>%
  group_by(clstr.name, gene) %>%
  # summarise(count.sum = sum(count.norm)) %>%
  summarise(count.sum = sum(count.raw)) %>%
  group_by(clstr.name) %>%
  mutate(count.norm = 10^6 * count.sum / sum(count.sum)) %>%
  rowwise() %>%
  mutate(logcount.norm = log2(count.norm + 1)) %>%
  group_by(gene) %>%
  mutate(zscore = scale(logcount.norm, center = TRUE, scale = TRUE))

count.mat.bulk <- count.long.merge %>%
  # dplyr::select(gene, clstr.name, logcount.norm) %>%
  # tidyr::spread(key = clstr.name, value = logcount.norm) %>%
  dplyr::select(gene, clstr.name, count.sum) %>%
  tidyr::spread(key = clstr.name, value = count.sum) %>%
  as.data.frame()


rownames(count.mat.bulk) <- count.mat.bulk$gene
count.mat.bulk$gene <- NULL

boxplot(count.mat.bulk, main = "Before quant norm")




# remove outliers, low number of genes?
cname.remove <- "neutro4"
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



rnames.orig <- rownames(count.mat.bulk)
cnames.orig <- colnames(count.mat.bulk)

count.mat.bulk <- preprocessCore::normalize.quantiles(as.matrix(count.mat.bulk), copy = TRUE)

rownames(count.mat.bulk) <- rnames.orig
colnames(count.mat.bulk) <- cnames.orig

boxplot(count.mat.bulk, main = "After quant norm")

dat.pca <- prcomp(t(count.mat.bulk), center = TRUE, scale. = TRUE)

dat.pca.proj <- as.data.frame(t(count.mat.bulk) %*% dat.pca$rotation %*% diag(dat.pca$sdev))
# dat.pca.proj <- as.data.frame(t(count.mat.bulk) %*% dat.pca$rotation)
colnames(dat.pca.proj) <- paste("PC", seq(ncol(dat.pca.proj)), sep = "")
dat.pca.proj$clstr.name <- rownames(dat.pca.proj)

ggplot(dat.pca.proj, aes(x = PC1, y = PC2, label = clstr.name)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_text_repel()


# Set up counts -----------------------------------------------------------

# remove bad genes
bad.genes <- "Gm1821"
count.mat.bulk <- count.mat.bulk[!grepl(bad.genes, rownames(count.mat.bulk)), ]
rownames(count.mat.bulk) <- sapply(rownames(count.mat.bulk), function(x) strsplit(x, "--")[[1]][[1]])

# make long
zscore.mat <- count.long.merge %>%
  dplyr::select(gene, clstr.name, zscore) %>%
  tidyr::spread(key = clstr.name, value = zscore) %>%
  as.data.frame()
rownames(zscore.mat) <- zscore.mat$gene
zscore.mat <- zscore.mat[!grepl(bad.genes, rownames(zscore.mat)), ]
zscore.mat$gene <- NULL
rownames(zscore.mat) <- sapply(rownames(zscore.mat), function(x) strsplit(x, "--")[[1]][[1]])

# dat.mat <- zscore.mat




# Do likelihoods ----------------------------------------------------------

# create likelihoods
jcutoff <- 1.8
plot(density(rowMeans(count.mat.bulk)))
abline(v = jcutoff)

dat.mat.filt <- count.mat.bulk[apply(count.mat.bulk, 1, max) > jcutoff, ]

genes.keep <- intersect(genes.all, rownames(dat.mat.filt))
terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])

dat.mat.filt <- dat.mat.filt[genes.keep, ]


# redefine genes that we keep
genes.filt <- intersect(rownames(dat.mat.filt), genes.keep)

dat.mat.filt <- dat.mat.filt[genes.filt, cname.keep.i]
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

# plot(probs.lst.filt$neutro5, probs.lst.filt$neutro3, pch = 20)
# text(probs.lst.filt$neutro5, probs.lst.filt$neutro3, labels = names(probs.lst.filt$neutro3))

# library(GGally)
# ggpairs(as.data.frame(dat.mat.filt))


# Prepare count mat -------------------------------------------------------

# count.filt <- count.dat$counts[terms.keep, ]
# count.filt <- count.dat$counts[terms.keep, ]
count.filt <- out.objs$count.mat[terms.keep, ]
rownames(count.filt) <- genes.keep

# Do fits -----------------------------------------------------------------

all.cells <- colnames(count.filt)
names(all.cells) <- all.cells

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
    # return(dmultinom(x = cell.vec, prob = jprob^(1/4), log = TRUE))
  })
})



# Output ------------------------------------------------------------------



# Summarize fits ----------------------------------------------------------

# calculate probability of model given data 
p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)

cell.counts <- Matrix::colSums(out.objs$count.mat) / 5
cell.counts.downsamp <- Matrix::colSums(out.objs$count.mat) / 5

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

# be stringent with the predictions 


LL.sum <- LL.dat %>%
  group_by(ctype.pred) %>%
  summarise(ncell = length(ctype.pred))

print(LL.sum)

p.filt <- log(0.99)
LL.dat <- LL.dat %>%
  rowwise() %>%
  mutate(ctype.stringent = ifelse(p.max >= p.filt, ctype.pred, NA))

LL.dat.merge <- left_join(dat.umap.pred.merged , LL.dat)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.umap.celltype <- ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = ctype.stringent)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey95") +
  theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.umap.celltype)

m.umap.celltype.nofilt <- ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey95") +
  theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.umap.celltype.nofilt)


ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.pred) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = LL.max/cell.size)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.stringent) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show enrihment before and after
m.enrichment.umap <- ggplot(subset(LL.dat.merge, !is.na(ctype.stringent)), aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~ctype.stringent)

print(m.enrichment.umap)

# Try integrating on scrna-seq spce ---------------------------------------

library(liger)

rna.clusts <- meta.neutro$clstr.name
names(rna.clusts) <- meta.neutro$cell

atac.clusts <- dat.umap.long.lda$louvain
names(atac.clusts) <- dat.umap.long.lda$cell

rna.mat <- count.long %>%
  dplyr::select(c(gene, cell, count.raw)) %>%
  tidyr::spread(key = cell, value = count.raw,) %>%
  as.data.frame()
rownames(rna.mat) <- rna.mat$gene
rna.mat$gene <- NULL

rna.mat <- rna.mat[!grepl(bad.genes, rownames(rna.mat)), ] 
rownames(rna.mat) <- sapply(rownames(rna.mat), function(x) strsplit(x, "--")[[1]][[1]])

bm.data = list(atac=count.filt, rna=rna.mat)
int.bm <- createLiger(bm.data)
int.bm <- normalize(int.bm)
int.bm <- selectGenes(int.bm, datasets.use = 2)
int.bm <- scaleNotCenter(int.bm)

sys.time(
  int.bm <- optimizeALS(int.bm, k=20)
)

int.pbmc <- runUMAP(int.bm, use.raw = T)
p1 <- plotByDatasetAndCluster(int.pbmc, return.plots = T)
print(p1[[1]])

int.bm <- quantileAlignSNF(int.bm)

int.bm <- runUMAP(int.bm, use.raw = T)
# p1 <- plotByDatasetAndCluster(int.bm, return.plots = T, clusters = rna.clusts)
p1 <- plotByDatasetAndCluster(int.bm, return.plots = T, clusters = atac.clusts)
print(p1[[1]])
print(p1[[2]])



# Get enrichment ----------------------------------------------------------

# calculate enrichment of each celltype
celltype.fc <- LL.dat.merge %>%
  filter(!is.na(ctype.stringent)) %>%
  group_by(is.stem, ctype.stringent) %>%
  summarise(cell.count = length(cell)) %>%
  group_by(is.stem) %>%
  mutate(cell.frac = cell.count / sum(cell.count))

jctypes <- unique(celltype.fc$ctype.stringent)
names(jctypes) <- jctypes
hyp.test.lst <- lapply(jctypes, DoFishersTestCtypeStringent, celltype.fc)

# plot odds ratio and p-value
hyp.test.dat <- lapply(jctypes, function(jctype){
  hyp.test <- hyp.test.lst[[jctype]]
  if (is.null(hyp.test)){
    warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
    return(NULL)
  }
  data.frame(OR = hyp.test$estimate, pval = hyp.test$p.value, ctype = jctype, stringsAsFactors = FALSE)
}) %>%
  bind_rows() %>%
  arrange(desc(OR)) %>%
  mutate(ctype = as.factor(ctype)) %>%
  dplyr::rename(ctype.stringent = ctype)

m.enrichment <- ggplot(hyp.test.dat, aes(x = log10(OR), y = -log10(pval), label = ctype.stringent)) + geom_point(size = 2.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text_repel(size = 4) +
  ggtitle("Fisher's exact test to quantify cell-type enrichment")
print(m.enrichment)

# sort by enrichment
celltype.fc.merge <- left_join(celltype.fc, hyp.test.dat)
m.barplot <- ggplot(celltype.fc.merge, aes(x = forcats::fct_reorder(.f = ctype.stringent, .x = OR, .desc = TRUE), y = cell.frac, group = is.stem, fill = is.stem)) + geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Celltype counts lineage-neg versus control, ordered by decreasing odds ratio")
print(m.barplot)





# 
# # Set up tables -----------------------------------------------------------
# 
# 
# # plot out
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 25
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# # umap.out <- umap(lsi.out$u, config = jsettings)
# 
# inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.Robj"
# load(inf.lda, v=T)
# 
# out.lda <- out.lda[[3]]
# 
# umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
# dat.umap.long.lda <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
# 
# dat.umap.long.lda$is.stem <- sapply(dat.umap.long.lda$cell, function(x) grepl("stem-cell", x))
# 
# ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2)) + geom_point() + 
#   facet_wrap(~is.stem) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# dat.umap.long.lda <- DoLouvain(topics.mat = posterior(out.lda)$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.lda)
# 
# ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
#   facet_wrap(~is.stem) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# inf.proj <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/projections/H3K4me1_matsMergedNonenriched.matsMergedEnriched.RData"
# load(inf.proj, v=T)
# 
# dat.pred <- predict(umap.out, out.lda.predict$topics)
# 
# dat.pred.long <- data.frame(cell = rownames(dat.pred), umap1 = dat.pred[, 1], umap2 = dat.pred[, 2], stringsAsFactors = FALSE) %>%
#   mutate(is.stem = TRUE)
# 
# dat.umap.pred.merged <- bind_rows(subset(dat.umap.long.lda, select = -louvain), dat.pred.long)
# 
# ggplot(dat.umap.pred.merged, aes(x = umap1, y = umap2)) + geom_point() + 
#   facet_wrap(~is.stem) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # Get raw counts ----------------------------------------------------------
# 
# inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedAll_2019-09-29.CountThres0.K-30_35_50.OutObjs.RData"
# load(inf.lda, v=T)
# 
# 
# 
# # Set up mean of clusters -------------------------------------------------
# 

# 
# 
# # Load bulk ---------------------------------------------------------------
# 
# inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
# dat <- fread(inf.bulkdat, sep = "\t")
# 
# colnames(dat) <- gsub(" ", "_", colnames(dat))
# dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
#   group_by(Gene_ID) %>%
#   mutate(FPKM = replace_na(FPKM, 0)) %>%
#   group_by(Gene_Name, CellType) %>%
#   summarise(FPKM = sum(FPKM)) %>%
#   rowwise() %>%
#   mutate(logFPKM = log2(FPKM + 1))
# 
# # normalize across samples?
# ggplot(dat.long, aes(x = CellType, y = logFPKM)) + geom_boxplot() +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# dat.mat <- tidyr::spread(dat.long %>%
#                            ungroup() %>%
#                            # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>%
#                            mutate(gene = Gene_Name) %>%
#                            dplyr::select(gene, CellType, logFPKM),
#                          key = CellType, value = logFPKM)  %>%
#   as.data.frame()
# rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL
# 
# cnames.tmp <- colnames(dat.mat)
# rnames.tmp <- rownames(dat.mat)
# dat.mat <- preprocessCore::normalize.quantiles(as.matrix(dat.mat), copy = TRUE)  # strong normalization,
# colnames(dat.mat) <- cnames.tmp
# rownames(dat.mat) <- rnames.tmp
# 
# boxplot(dat.mat)
# 
# dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "exprs", -gene) %>%
#   group_by(gene) %>%
#   mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))
# 
# dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>%
#   as.data.frame()
# rownames(dat.norm.zscore.mat) <- dat.norm.zscore.mat$gene
# dat.norm.zscore.mat$gene <- NULL
# 
# annot.out <- AnnotateBins(terms.mat = out.objs$tm.result$terms)
# 
# terms.sum <- annot.out$terms.filt %>%
#   group_by(gene) %>%
#   dplyr::filter(rnk == min(rnk))
# 
# term2gene <- hash(terms.sum$term, terms.sum$gene)
# gene2term <- hash(terms.sum$gene, terms.sum$term)
# 
# genes.keep <- rownames(dat.mat)
# 
# top.nterms <- 150
# terms.all <- unique(subset(terms.sum %>% arrange(desc(weight)), rnk <= top.nterms)$term)
# genes.all <- unlist(sapply(terms.all, function(x) term2gene[[x]]))
# 
# # Do likelihoods ----------------------------------------------------------
# 
# # create likelihoods
# jcutoff <- 1.8
# plot(density(rowMeans(dat.mat)))
# abline(v = jcutoff)
# 
# dat.mat.filt <- dat.mat[apply(dat.mat, 1, max) > jcutoff, ]
# 
# genes.keep <- intersect(genes.all, rownames(dat.mat.filt))
# terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])
# 
# dat.mat.filt <- dat.mat.filt[genes.keep, ]
# 
# 
# 
# # redefine genes that we keep
# genes.filt <- intersect(rownames(dat.mat.filt), genes.keep)
# 
# # dat.mat.filt <- dat.mat.filt[genes.filt, ]
# # try using zscore
# dat.mat.filt <- 2^dat.norm.zscore.mat[genes.filt, ]  # exp is worse than 2
# 
# # dat.mat.filt <- dat.norm.zscore.mat[genes.filt, ]
# # renormalize so no zeros??
# # dat.mat.filt <- sweep(dat.mat.filt, MARGIN = 2, STATS = apply(dat.mat.filt, 2, min), FUN = "-")
# # dat.mat.filt <- as.matrix(dat.mat.filt)
# 
# 
# # handle zeros
# zero.fill <- min(as.matrix(dat.mat.filt)[which(as.matrix(dat.mat.filt) != 0)])
# dat.mat.filt[which(dat.mat.filt == 0)] <- zero.fill
# 
# # make likelihoods
# probs.lst.filt <- as.list(as.data.frame(dat.mat.filt))
# # name the list just to be safe
# probs.lst.filt <- lapply(probs.lst.filt, function(x){
#   names(x) <- rownames(dat.mat.filt)
#   return(x)
# }) 
# 
# 
# # Prepare count mat -------------------------------------------------------
# 
# count.filt <- count.dat$counts[terms.keep, ]
# rownames(count.filt) <- genes.keep
# 
# # Do fits -----------------------------------------------------------------
# 
# all.cells <- colnames(count.dat$counts)
# names(all.cells) <- all.cells
# 
# LL.ctype.lst <- lapply(all.cells, function(cell.name){
#   cell.vec <- count.filt[genes.filt, cell.name]
#   # cell.vec <- cell.vec[which(cell.vec > 0)]
#   LL.vec <- sapply(probs.lst.filt, function(jprob){
#     # jprob <- jprob[names(cell.vec)]
#     assertthat::assert_that(all(names(cell.vec) == names(jprob)))
#     # return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
#     # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
#     # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
#     # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
#     return(dmultinom(x = cell.vec, prob = jprob^(1/4), log = TRUE))
#   })
# })
# 
# 
# 
# # Output ------------------------------------------------------------------
# 
# 
# 
# # Summarize fits ----------------------------------------------------------
# 
# # calculate probability of model given data 
# p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)
# 
# cell.counts <- Matrix::colSums(count.dat$counts) / 5
# cell.counts.downsamp <- Matrix::colSums(count.dat$counts) / 5
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
# ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = ctype.stringent)) + geom_point() + 
#   scale_color_manual(values = cbPalette) +
#   theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
#   scale_color_viridis_c() + 
#   facet_wrap(~ctype.stringent) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = LL.max/cell.size)) + geom_point() + 
#   scale_color_viridis_c() + 
#   facet_wrap(~ctype.stringent) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # show enrihment before and after
# ggplot(subset(LL.dat.merge, !is.na(ctype.stringent)), aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + 
#   scale_color_manual(values = cbPalette) +
#   theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~ctype.stringent)
# 
# 
# 
# # Get enrichment ----------------------------------------------------------
# 
# # calculate enrichment of each celltype
# celltype.fc <- LL.dat.merge %>%
#   filter(!is.na(ctype.stringent)) %>%
#   group_by(is.stem, ctype.stringent) %>%
#   summarise(cell.count = length(cell)) %>%
#   group_by(is.stem) %>%
#   mutate(cell.frac = cell.count / sum(cell.count))
# 
# jctypes <- unique(celltype.fc$ctype.stringent)
# names(jctypes) <- jctypes
# hyp.test.lst <- lapply(jctypes, DoFishersTestCtypeStringent, celltype.fc)
# 
# # plot odds ratio and p-value
# hyp.test.dat <- lapply(jctypes, function(jctype){
#   hyp.test <- hyp.test.lst[[jctype]]
#   if (is.null(hyp.test)){
#     warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
#     return(NULL)
#   }
#   data.frame(OR = hyp.test$estimate, pval = hyp.test$p.value, ctype = jctype, stringsAsFactors = FALSE)
# }) %>%
#   bind_rows() %>%
#   arrange(desc(OR)) %>%
#   mutate(ctype = as.factor(ctype)) %>%
#   dplyr::rename(ctype.stringent = ctype)
# 
# m.enrichment <- ggplot(hyp.test.dat, aes(x = log10(OR), y = -log10(pval), label = ctype.stringent)) + geom_point(size = 2.5) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   geom_text_repel(size = 4) +
#   ggtitle("Fisher's exact test to quantify cell-type enrichment")
# print(m.enrichment)
# 
# # sort by enrichment
# celltype.fc.merge <- left_join(celltype.fc, hyp.test.dat)
# m.barplot <- ggplot(celltype.fc.merge, aes(x = forcats::fct_reorder(.f = ctype.stringent, .x = OR, .desc = TRUE), y = cell.frac, group = is.stem, fill = is.stem)) + geom_bar(stat = "identity", position = "dodge") +
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   ggtitle("Celltype counts lineage-neg versus control, ordered by decreasing odds ratio")
# print(m.barplot)
# 
# 
# 



# 
# # Why monocytes are highlighting two cells ?  -----------------------------
# 
# cells.mono.up <- subset(LL.dat.merge, ctype.stringent == "monocyte" & umap1 < 0 & umap2 > -2.5)$cell
# cells.mono.down <- subset(LL.dat.merge, ctype.stringent == "monocyte" & umap1 < 0 & umap2 < -2.5)$cell
# 
# 
# 
# counts.up <- Matrix::rowSums(count.filt[genes.filt, cells.mono.up])
# counts.down <- Matrix::rowSums(count.filt[genes.filt, cells.mono.down])
# 
# counts.mono <- data.frame(gene = names(counts.up), up = counts.up, down = counts.down, stringsAsFactors = FALSE)
# 
# ggplot(counts.mono, aes(x = up, y = down, label = gene)) + 
#   geom_point() + geom_text() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_abline(slope = 1, intercept = 0)
# 
# 
# counts.mono <- counts.mono %>%
#   rowwise() %>%
#   mutate(FC.down.vs.up = log2(down / up), 
#          mean.counts = mean(c(up, down))) %>%
#   mutate(gene.lab = ifelse(abs(FC.down.vs.up) > 1.25, gene, NA))
# 
# ggplot(counts.mono, aes(x = mean.counts, y = FC.down.vs.up, label = gene.lab)) + 
#   geom_text_repel() + 
#   geom_point() + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_hline(yintercept = 0)
# 
# 
# 
# # Test a cell -------------------------------------------------------------
# 
# plot(probs.lst.filt$granulocyte, probs.lst.filt$monocyte, pch = 20)
# text(probs.lst.filt$granulocyte, probs.lst.filt$monocyte, labels = names(probs.lst.filt[[1]]))
# 
# # jcell <- cells.mono.up[[1]]
# jcell <- cells.mono.down[[1]]
# # jcell <- (subset(LL.dat.merge, ctype.stringent == "granulocyte") %>% filter(umap2 == min(umap2)))$cell
# 
# # compare likelihoods
# Lvec <- sort(LL.ctype.lst[[jcell]], decreasing = TRUE)
# Pvec <- sort(p.ctype.lst[[jcell]], decreasing = TRUE)
# 
# plot(x = seq(length(Lvec)), y = Lvec)
# text(x = seq(length(Lvec)), y = Lvec, labels = names(Lvec))
# 
# # plot cell on umap
# ggplot(LL.dat.merge %>% mutate(is.cell = cell %in% jcell), aes(x = umap1, y = umap2, color = is.cell, size = is.cell)) + geom_point() +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # look at raw vector see why
# cell.vec <- count.filt[genes.filt, jcell][names(probs.lst.filt[[1]])]
# 
# # jgenes.remove <- c("")
# jgenes.remove <- names(which(cell.vec == 0))
# par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# jctypes <- c("monocyte", "granulocyte", "megakaryocyte", "nucleate_erythrocyte")
# # jctype <- "monocyte"
# for (jctype in jctypes){
#   print(jctype)
#   # ptmp <- probs.lst.filt[[jctype]] ^ (1 / 4)
#   ptmp <- log2(probs.lst.filt[[jctype]] + 1)
#   # ptmp <- sqrt(probs.lst.filt[[jctype]])
#   # ptmp <- probs.lst.filt[[jctype]]
#   ptmp <- ptmp / sum(ptmp)
#   jgenes.filt <- names(ptmp)[which(!names(ptmp) %in% jgenes.remove)]
#   
#   # L <- dmultinom(x = cell.vec[jgenes.filt], prob = ptmp[jgenes.filt], log = TRUE)
#   L <- dmultinom(x = cell.vec, prob = ptmp, log = TRUE)
#   # plot(x = cell.vec, y = ptmp, main = paste(jctype, signif(L, digits = 4)), pch = 20, xlab = jcell, ylab = jctype)
#   # text(x = cell.vec, y = ptmp, labels = names(cell.vec))
#   plot(y = cell.vec, x = ptmp, main = paste(jctype, signif(L, digits = 4)), pch = 20, xlab = jcell, ylab = jctype)
#   text(y = cell.vec, x = ptmp, labels = names(cell.vec))
#   # plot(y = cell.vec[jgenes.filt], x = ptmp[jgenes.filt], main = paste(jctype, signif(L, digits = 4)), pch = 20, xlab = jcell, ylab = jctype)
#   # text(y = cell.vec[jgenes.filt], x = ptmp[jgenes.filt], labels = names(cell.vec))
# }
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# 
# 
# # what is exprs of top granu genes across counts?
# jsub <- count.filt[genes.filt, ]
