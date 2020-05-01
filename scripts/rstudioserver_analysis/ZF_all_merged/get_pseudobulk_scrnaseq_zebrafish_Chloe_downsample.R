# Jake Yeung
# Date of Creation: 2020-04-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/get_pseudobulk_scrnaseq_zebrafish_Chloe_downsample.R
# Filter out bad cells? do pseudobulk but also downsample 


rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(xlsx)
library(Seurat)
library(hash)

library(DropletUtils)

library(scchicFuncs)
library(JFuncs)

Vectorize(plog2p <- function(p){
  return(ifelse(p == 0, 0, p * log2(p)))
}, vectorize.args = "p")

CalculateEntropy <- function(p, normalize.p = FALSE){
  if (normalize.p){
    p <- p / sum(p)
  }
  S <- -sum(plog2p(p))
  return(S)
}


# Load marker genes  ------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/For_Jake"


jlogfc.threshold <- 0
jmin.pct <- 0
jtest <- "poisson"
# jtest <- "DESeq2"
# jtest <- "MAST"
# jtest <- "negbinom"

outdir <- paste0("/home/jyeung/data/from_rstudioserver/zebrafish.", jtest, ".", Sys.Date())  # run over night
dir.create(outdir)
outrds <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.rds"))
outrds.ctypefilt <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.ctypefilt.rds"))
outpdf <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.pdf"))

inf <- file.path(indir, paste0("the_massive_complete_zf_dataset.csv.gz"))
inf.meta <- file.path(indir, "tsne_clusterID_zebrafish_GateID_dataset.csv")
inf.meta2 <- file.path(indir, "cluster_info_JY_edited.xlsx")


dat <- as.data.frame(fread(inf, stringsAsFactors = FALSE))
colnames(dat)[[1]] <- "gene"
rownames(dat) <- dat$gene
dat$gene <- NULL

mat <- dat

meta <- fread(inf.meta, stringsAsFactors = TRUE)
meta2 <- read.xlsx2(inf.meta2, sheetIndex = 1, header = FALSE, colClasses = c("integer", "character")); colnames(meta2) <- c("ClusterID", "celltype")

colnames(meta) <- c("rname", "V1", "V2", "ClusterID", "experi")

meta <- left_join(meta, meta2)

meta$clusterplate <- as.character(paste(meta$celltype, meta$experi, sep = "_"))
meta$celltype.merge <- sapply(as.character(meta$celltype), function(x) ifelse(x %in% c("monocytes", "neutrophils"), "granulocytes", x))

table(meta$clusterplate)

# Check depth per cell ----------------------------------------------------


pdf(file = outpdf, useDingbats = FALSE)

mat.filt <- mat[, meta$rname]

plot(density(colSums(mat.filt)), log = "x")

ngenes <- apply(mat.filt, MARGIN = 2, function(jcol) Matrix::nnzero(jcol))
nreads <- apply(mat.filt, MARGIN = 2, function(jcol) sum(jcol))

plot(density(ngenes), log = "x")
plot(density(nreads), log = "x")


# Filter out bad cells ----------------------------------------------------

min.ngenes <- 100
min.nreads <- 300

cells.keep.i <- which(ngenes >= min.ngenes & nreads >= min.nreads)
cells.keep <- names(ngenes)[cells.keep.i]

print(length(cells.keep))

# Do pseudobulks on cleaner cells -----------------------------------------

mat.filt2 <- mat.filt[, cells.keep]

cnames.keep.lst <- lapply(split(meta, meta$celltype), function(x){
  as.character(x$rname)
})

cnames.keep.lst.ctypefilt <- lapply(split(meta, meta$celltype.merge), function(x){
  as.character(x$rname)
})

mat.pbulk <- SumAcrossClusters(mat.filt2, cnames.keep.lst)
mat.pbulk.ctypefilt <- SumAcrossClusters(mat.filt2, cnames.keep.lst.ctypefilt)

mat.pbulk <- do.call(cbind, mat.pbulk)
mat.pbulk.ctypefilt <- do.call(cbind, mat.pbulk.ctypefilt)

boxplot(log2(mat.pbulk))
boxplot(log2(mat.pbulk.ctypefilt))

# downsample

set.seed(0)
jprop <- min(colSums(mat.pbulk)) / colSums(mat.pbulk)
jprop.ctypefilt <- min(colSums(mat.pbulk.ctypefilt)) / colSums(mat.pbulk.ctypefilt)
mat.pbulk.ds <- DropletUtils::downsampleMatrix(mat.pbulk, jprop, bycol=TRUE)
mat.pbulk.ds.ctypefilt <- DropletUtils::downsampleMatrix(mat.pbulk.ctypefilt, jprop.ctypefilt, bycol=TRUE)

# remove genes with zero variance
# pbulk.long <- data.frame(gene = rownames(mat.pbulk.ds.filt), mat.pbulk.ds.filt, stringsAsFactors = TRUE) %>%
#   reshape2::melt(., id.vars = "gene", variable.name = "pbulk", value.name = "counts") %>% 
#   group_by(gene) %>%
#   mutate(log2p1counts = log2(counts + 1), 
#          log2fc = log2p1counts - mean(log2p1counts), 
#          log2zscore = log2fc / sd(log2p1counts))

jvar.min <- 0
jmean.min <- 4

plot(density(log2(mat.pbulk.ds + 1))); abline(v = jmean.min)
plot(density(log2(mat.pbulk.ds.ctypefilt + 1))); abline(v = jmean.min)

mat.pbulk.ds.filt <- mat.pbulk.ds[which(apply(mat.pbulk.ds, MARGIN = 1, function(jrow) var(jrow)) > jvar.min), ]
mat.pbulk.ds.filt <- mat.pbulk.ds.filt[which(apply(mat.pbulk.ds.filt, MARGIN = 1, function(jrow) mean(jrow)) > jmean.min), ]
mat.pbulk.ds.ctypefilt.filt <- mat.pbulk.ds.ctypefilt[which(apply(mat.pbulk.ds.ctypefilt, MARGIN = 1, function(jrow) var(jrow)) > jvar.min), ]
mat.pbulk.ds.ctypefilt.filt <- mat.pbulk.ds.ctypefilt.filt[which(apply(mat.pbulk.ds.ctypefilt.filt, MARGIN = 1, function(jrow) mean(jrow)) > jmean.min), ]

plot(density(log2(mat.pbulk.ds + 1)))
plot(density(log2(mat.pbulk.ds.filt + 1)))
plot(density(log2(mat.pbulk.ds.ctypefilt.filt + 1)))

pca.out <- prcomp(t(log2(mat.pbulk.ds.filt + 1)), center = TRUE, scale. = TRUE)

dat.pca <- data.frame(cell = rownames(pca.out$x), pca.out$x, stringsAsFactors = FALSE)

ggplot(dat.pca, aes(x = PC1, y = PC2, label = cell)) + geom_point()  + geom_text() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

ggplot(dat.pca, aes(x = PC2, y = PC3, label = cell)) + geom_point()  + geom_text() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

ggplot(dat.pca, aes(x = PC3, y = PC4, label = cell)) + geom_point()  + geom_text() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)


# Estimate log2FC and zscores ---------------------------------------------

pbulk.long <- data.frame(gene = rownames(mat.pbulk.ds.filt), mat.pbulk.ds.filt, stringsAsFactors = TRUE) %>%
  reshape2::melt(., id.vars = "gene", variable.name = "pbulk", value.name = "counts") %>% 
  group_by(gene) %>%
  mutate(log2p1counts = log2(counts + 1), 
         log2fc = log2p1counts - mean(log2p1counts), 
         zscore = log2fc / sd(log2p1counts))


pbulk.ctypefilt.long <- data.frame(gene = rownames(mat.pbulk.ds.ctypefilt.filt), mat.pbulk.ds.ctypefilt.filt, stringsAsFactors = TRUE) %>%
  reshape2::melt(., id.vars = "gene", variable.name = "pbulk", value.name = "counts") %>% 
  group_by(gene) %>%
  mutate(log2p1counts = log2(counts + 1), 
         log2fc = log2p1counts - mean(log2p1counts), 
         zscore = log2fc / sd(log2p1counts))


 
# Write outputs -----------------------------------------------------------

outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_pseudobulk_scrnaseq_downsampled.", Sys.Date(), ".RData")
save(mat.pbulk, mat.pbulk.ctypefilt, mat.pbulk.ds, mat.pbulk.ds.ctypefilt, pbulk.long, pbulk.ctypefilt.long, file = outf)



# Run seurat? -------------------------------------------------------------


jstart <- Sys.time()


zf <- CreateSeuratObject(counts = mat.filt2, project = "zf", min.cells = 10, min.features = 400)

plot2 <- FeatureScatter(zf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot2)


zf <- FindVariableFeatures(zf, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(zf), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(zf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
CombinePlots(plots = list(plot1, plot2))

# add labels?
zf@meta.data$cell <- rownames(zf@meta.data)
jhash <- hash(meta$rname, meta$experi)
jhash2 <- hash(meta$rname, meta$ClusterID)
jhash3 <- hash(meta$rname, as.character(meta$celltype))
jhash4 <- hash(meta$rname, as.character(meta$celltype.merge))
zf@meta.data$experi <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash[[x]], NA)
})
zf@meta.data$clusterid <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash2[[x]]), jhash2[[x]], NA)
})
zf@meta.data$celltype <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash3[[x]], NA)
})
zf@meta.data$celltype.merge <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash4[[x]], NA)
})

all.genes <- rownames(zf)
zf <- ScaleData(zf, features = all.genes)
zf <- RunPCA(zf, features = VariableFeatures(object = zf))
zf <- RunUMAP(zf, dims = 1:10)
DimPlot(zf, reduction = "umap", label = FALSE, group.by = "clusterid", na.value = "#C0C0C0") + scale_color_viridis_d()

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
DimPlot(zf, reduction = "umap", label = FALSE, group.by = "celltype", na.value = "#C0C0C0", cols = cbPalette) + ggtitle("UMAP from log and scale the data")

dev.off()

# Do DE expression  -------------------------------------------------------


# Get markers -------------------------------------------------------------


Idents(zf) <- zf$celltype.merge
print("Running find markers ctypefilt")
de.output.ctypefilt.lst <- FindAllMarkers(zf, only.pos = FALSE, logfc.threshold = jlogfc.threshold, min.pct = jmin.pct, test.use = jtest)
print("done Running find markers")
saveRDS(de.output.ctypefilt.lst, file = outrds.ctypefilt)
print(Sys.time() - jstart)




Idents(zf) <- zf$celltype
print("Running find markers")
de.output.lst <- FindAllMarkers(zf, only.pos = FALSE, logfc.threshold = jlogfc.threshold, min.pct = jmin.pct, test.use = jtest)
print("done Running find markers")
saveRDS(de.output.lst, file = outrds)
print(Sys.time() - jstart)
