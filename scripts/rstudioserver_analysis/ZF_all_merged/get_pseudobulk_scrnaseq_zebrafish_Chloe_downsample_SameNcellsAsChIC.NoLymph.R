# Jake Yeung
# Date of Creation: 2020-04-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/get_pseudobulk_scrnaseq_zebrafish_Chloe_downsample_SameNcellsAsChIC.R
# 



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


# Load chic reference -----------------------------------------------------


# keep just a few ctypes?
# ctypes.keep <- c("granulocytes", "erythrocytes", "HSPCs", "lymphocytes")
ctypes.keep <- c("granulocytes", "erythrocytes", "HSPCs")



# load K4me1 annotations from WKM
jmark.ref <- "H3K4me1"
inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark.ref, ".keepn_150.final.ClusterTables.txt")
assertthat::assert_that(file.exists(inf.annot.louv))
inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark.ref, ".2020-04-13.txt")
assertthat::assert_that(file.exists(inf.annot.glmpca))

annot.louv <- fread(inf.annot.louv)
annot.glmpca <- fread(inf.annot.glmpca)

annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")

annot.glmpca <- fread(inf.annot.glmpca)
annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
  rowwise() %>%
  mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
  mutate(cluster = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
  ungroup() %>%
  filter(!cluster %in% c("Unknown", "thrombo"))  # remove the small cluster Unknown

# count number of cells per cluster, remoing lymph?
jsplit <- lapply(split(annot.glmpca.filt, f = annot.glmpca.filt$cluster), function(x){
  return(subset(x, select = c(cell, cluster)))
}) %>%
  bind_rows()
jsplit.sum <- jsplit %>%
  group_by(cluster) %>%
  summarise(ncells = length(cell)) %>%
  ungroup() %>%
  mutate(nfrac = ncells / sum(ncells),
         prop = min(nfrac) / nfrac,
         cluster = gsub("eryth", "erythrocytes", cluster),
         cluster = gsub("HSC", "HSPCs", cluster),
         cluster = gsub("monocyte", "granulocytes", cluster),
         cluster = gsub("lymph", "lymphocytes", cluster))

print(jsplit.sum$cluster)


# Load marker genes  ------------------------------------------------------

make.plots <- FALSE

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/For_Jake"


jlogfc.threshold <- 0
jmin.pct <- 0
jtest <- "poisson"
# jtest <- "DESeq2"
# jtest <- "MAST"
# jtest <- "negbinom"

outdir <- paste0("/home/jyeung/data/from_rstudioserver/zebrafish.", jtest, ".SameNbrCells2.NoLymph.", Sys.Date())  # run over night
dir.create(outdir)

outf <- file.path(outdir, paste0("WKM_pseudobulk_scrnaseq_downsampled.", Sys.Date(), ".SameNbrCells.NoLymph.RData"))
outrds <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.SameNbrCells.rds"))
outrds.ctypefilt <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.ctypefilt.SameNbrCells.NoLymph.rds"))
outpdf <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.SameNbrCells.pdf"))

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

if (make.plots){
  pdf(file = outpdf, useDingbats = FALSE)
}

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
cells.keep <- names(ngenes)[cells.keep.i]  # need to define again later after downsampling





print(length(cells.keep))

# Do pseudobulks on cleaner cells -----------------------------------------

mat.filt2 <- mat.filt[, cells.keep]

# redefine meta
meta <- subset(meta, rname %in% cells.keep)

cnames.keep.lst <- lapply(split(meta, meta$celltype), function(x){
  as.character(x$rname)
})

lapply(cnames.keep.lst, length)

cnames.keep.lst.ctypefilt <- lapply(split(meta, meta$celltype.merge), function(x){
  as.character(x$rname)
})

ncells.unfilt <- sapply(cnames.keep.lst, length)

# downsample simply to lowest
N.lowest <- min(ncells.unfilt)
cnames.keep.lst <- lapply(cnames.keep.lst, function(xsub){
  sample(xsub, size = N.lowest, replace = FALSE)
})
lapply(cnames.keep.lst, length)



for (jname in names(cnames.keep.lst.ctypefilt)){
  if (!jname %in% ctypes.keep){
    cnames.keep.lst.ctypefilt[[jname]] <- NULL
  }
}

ncells.scrnaseq <- sapply(cnames.keep.lst.ctypefilt, length) 

# downsample so that same proportion of celltypes as in WKM chic data? or equal? 
# target proportions: 
prop.target <- jsplit.sum
jprop.min <- subset(prop.target, cluster == "erythrocytes")$nfrac
props.after <- prop.target$nfrac; names(props.after) <- prop.target$cluster
total.cells.after <- min(ncells.scrnaseq) / jprop.min
ncells.after <- round(total.cells.after * props.after)

# downsample cells to target 
jctypes <- names(cnames.keep.lst.ctypefilt)
names(jctypes) <- jctypes

print("Downsampling these celltypes:")
print(jctypes)

cnames.keep.lst.ctypefilt.downsamp <- lapply(jctypes, function(jctype){
  # ncells after
  print(jctype)
  xsamp <- sample(x = cnames.keep.lst.ctypefilt[[jctype]], size = ncells.after[[jctype]])
  return(xsamp)
})

N <- length(unlist(cnames.keep.lst.ctypefilt.downsamp))
(lapply(cnames.keep.lst.ctypefilt.downsamp, function(x) length(x) / N))
print(prop.target)

# make hat the new ctpyefilt
cnames.keep.lst.ctypefilt <- cnames.keep.lst.ctypefilt.downsamp


cells.keep.ds <- unlist(cnames.keep.lst.ctypefilt)

print(cells.keep.ds)



# Continue as normal ------------------------------------------------------


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

save(mat.pbulk, mat.pbulk.ctypefilt, mat.pbulk.ds, mat.pbulk.ds.ctypefilt, pbulk.long, pbulk.ctypefilt.long, file = outf)



# Run seurat? -------------------------------------------------------------


jstart <- Sys.time()

# use only cells.keep.ds....
# replace mat.fit2 and continue as normal 
mat.filt2 <- mat.filt2[, cells.keep.ds]

print("Creating seurat object of size:")
print(dim(mat.filt2))

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

if (make.plots){
  dev.off()
}

# Do DE expression  -------------------------------------------------------

# CalculateVarRaw



# Get markers -------------------------------------------------------------


Idents(zf) <- zf$celltype
print("Running find markers")
de.output.lst <- FindAllMarkers(zf, only.pos = FALSE, logfc.threshold = jlogfc.threshold, min.pct = jmin.pct, test.use = jtest)
print("done Running find markers")
saveRDS(de.output.lst, file = outrds)
print(Sys.time() - jstart)



Idents(zf) <- zf$celltype.merge
print("Running find markers ctypefilt")
de.output.ctypefilt.lst <- FindAllMarkers(zf, only.pos = FALSE, logfc.threshold = jlogfc.threshold, min.pct = jmin.pct, test.use = jtest)
print("done Running find markers")
saveRDS(de.output.ctypefilt.lst, file = outrds.ctypefilt)
print(Sys.time() - jstart)


