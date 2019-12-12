# Jake Yeung
# Date of Creation: 2019-12-06
# File: ~/projects/scchic/scripts/scripts_analysis/public_data_transcriptome_grun/get_pseudobulk_from_grun.R
# Pseudobulk from Grun


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

jgene <- grep(pattern = "Elane", rownames(dat.filt), value = TRUE)

FeaturePlot(bm, features = jgene, pt.size = 0.2, ncol = 3)

# meta <- left_join(meta, data.frame(cell = rownames(bm@reductions$umap@cell.embeddings), bm@reductions$umap@cell.embeddings, stringsAsFactors = FALSE))

# Downstream --------------------------------------------------------------

count.long <- as.data.frame(bm@assays$RNA@counts) %>%
  ungroup() %>%
  mutate(gene = rownames(as.data.frame(bm@assays$RNA@counts))) %>%
  tidyr::gather(key = "cell", value = "count.raw", -gene)



# rename clusters
# clstrs.keep <- c(1, 8, 9, 5, 11, 28, 7, 35, 3, 6, 32, 4, 2, 14, 22, 13, 12, 10, 29)
# clstrs.rename <- c("HSC", "Prog8", "")
clstrs.keep <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 22, 28, 29, 32, 35)
clstrs.rename <- c("HSC", "Neutrophils", "NeutrophilsProg", "Megakaryocytes", "Macrophages", 
                   "Erythroblasts", "Erythroblasts", "Erythroblasts", "ErythroblastsProg", "Eosinophils", 
                   "NeutrophilsProg2", "Neutrophils", "Erythroblasts", "Neutrophils", "Neutrophils",
                   "Macrophages", "Eosinophils", "PhasgoEry", "Bcells")

clstrs.dat <- data.frame(cluster = clstrs.keep, clstr.name = clstrs.rename, stringsAsFactors = FALSE)

# clstrs.keep <- c(1, 11, 3, 2, 14, 12)
# clstrs.rename <- c("stem", "neutro1", "neutro2", "neutro3", "neutro4", "neutro5")
# clstrs.hash <- hash(clstrs.keep, clstrs.rename)
# clstrs.hash[["Other"]] <- "Other"

# average across  clusters
meta.ctype <- left_join(meta, clstrs.dat) %>%
  dplyr::filter(!is.na(clstr.name))

rownames(meta.ctype) <- meta.ctype$cell


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


# sum across cells
count.long.merge <- left_join(count.long, meta.ctype) %>%
  group_by(clstr.name, gene) %>%
  summarise(count.sum = sum(count.raw)) %>%
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


# # remove outliers, low number of genes?
# cname.remove <- "neutro4"
# cname.keep.i <- which(!colnames(count.mat.bulk) %in% cname.remove)
# count.mat.bulk <- count.mat.bulk[, cname.keep.i]


# Remove lowly expressed genes --------------------------------------------

# plot(density(log2(unlist(subset(count.mat.bulk, select = -gene)) + 1)))

count.means <- count.long.merge %>%
  group_by(gene) %>%
  summarise(count.norm.med = quantile(logcount.norm, 0.5))

jcutoff <- 3
ggplot(count.means, aes(x = count.norm.med)) + geom_density() + 
  geom_vline(xintercept = jcutoff) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

genes.filt <- subset(count.means, count.norm.med > jcutoff)$gene

# count.mat.bulk.mat <- subset(count.mat.bulk, select = -gene)
# count.means <- rowSums(count.mat.bulk)
# count.mat.bulk[1:5, 1:5]


# Do DESeq2 ---------------------------------------------------------------


# do DESeq2 variance stabilization 

metadata <- data.frame(ctype = colnames(count.mat.bulk), stringsAsFactors = FALSE)
rownames(metadata) <- colnames(count.mat.bulk)

count.mat.bulk.filt <- count.mat.bulk[genes.filt, ]
count.mat.bulk.int <- matrix(as.integer(as.matrix(count.mat.bulk.filt)), 
                             nrow = nrow(count.mat.bulk.filt), ncol = ncol(count.mat.bulk.filt), 
                             dimnames = list(rownames(count.mat.bulk.filt), colnames(count.mat.bulk.filt)))
dds <- DESeqDataSetFromMatrix(countData = count.mat.bulk.int, colData = metadata, design = ~1)
vsd <- vst(dds)

boxplot(assay(vsd))



rnames.orig <- rownames(count.mat.bulk.int)
cnames.orig <- colnames(count.mat.bulk.int)

count.mat.bulk.qn <- preprocessCore::normalize.quantiles(assay(vsd), copy = TRUE)

rownames(count.mat.bulk.qn) <- rnames.orig
colnames(count.mat.bulk.qn) <- cnames.orig

boxplot(count.mat.bulk.qn, main = "After quant norm")

dat.pca <- prcomp(t(count.mat.bulk.qn), center = TRUE, scale. = TRUE)

dat.pca.proj <- as.data.frame(t(count.mat.bulk.qn) %*% dat.pca$rotation)
# dat.pca.proj <- as.data.frame(t(count.mat.bulk) %*% dat.pca$rotation)
colnames(dat.pca.proj) <- paste("PC", seq(ncol(dat.pca.proj)), sep = "")
dat.pca.proj$clstr.name <- rownames(dat.pca.proj)

ggplot(dat.pca.proj, aes(x = PC1, y = PC2, label = clstr.name)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_text_repel()


# Save pseudobulks to output  ---------------------------------------------

save(count.mat.bulk.qn, count.long.merge, bm, file = "/Users/yeung/data/scchic/public_data/bone_marrow_grun/pseudobulk_celltypes_qnorm.RData")



# Entropy of pseudobulks?  ------------------------------------------------

library(DropletU)

bm <- NormalizeData(bm, normalization.method = "LogNormalize", scale.factor = 10000)
bm <- FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(bm)
bm <- ScaleData(bm, features = all.genes)

# S.vec <- apply(bm@assays$SCT@scale.data, 2, function(jcell) CalculateEntropy(2^jcell, normalize.p = TRUE))
# S.vec <- apply(bm@assays$RNA@data, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
S.vec <- apply(bm@assays$RNA@counts, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
jmeta <- data.frame(S = S.vec)
rownames(jmeta) <- names(S.vec)

bm <- AddMetaData(object = bm, metadata = jmeta, col.name = "entropy")

FeaturePlot(bm, features = 'entropy', pt.size = 3) + scale_color_viridis_c()

DimPlot(bm, label = TRUE, group.by = "cluster", cols = cbPalette)

meta.ctype.full <- as.data.frame(subset(left_join(meta, clstrs.dat), select = c(clstr.name, cell)))
rownames(meta.ctype.full) <- meta.ctype.full$cell
meta.ctype.full$cell <- NULL

bm <- AddMetaData(bm, metadata = meta.ctype.full, col.name = "clstr.name")

DimPlot(bm, label = TRUE, group.by = "clstr.name", cols = cbPalette, pt.size = 3) + ggtitle("Gruen et al. 2016",)



