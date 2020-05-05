# Jake Yeung
# Date of Creation: 2020-05-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/2-explore_Giladi_et_al_write_gene_tables_from_celltypefilt.KeepFourCelltypes.Downsample.R
# Try to get same number of cells in each cluster, use poisson, and keep only Bcells, erythroblasts, HSPCs, and "granulocytes"


rm(list=ls())

jstart <- Sys.time()

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(xlsx)
library(JFuncs)
library(scchicFuncs)

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



  
# Functions ---------------------------------------------------------------


ReadGiladi <- function(inf, remove.first.col = TRUE){
  dat <- fread(inf) %>% 
    dplyr::rename(gene = V1)
  if (remove.first.col){
    dat$gene <- NULL
  }
  return(dat)
}


# Load marker genes  ------------------------------------------------------

make.plots <- TRUE
jmeth <- "poisson"
WithVpreb <- FALSE
outrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.meth_", jmeth, ".Downsampled.FourCtypes.WithVpreb_", WithVpreb, ".rds")
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt..meth_", jmeth, ".Downsampled.FourCtypes.WithVpreb_", WithVpreb, ".pdf")


if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

inf.markers <- "/home/jyeung/data/from_cluster/public_data/Giladi_et_al_2018/41556_2018_121_MOESM4_ESM.markergenes.xlsx"
assertthat::assert_that(file.exists(inf.markers))


markers <- xlsx::read.xlsx(inf.markers, sheetIndex = 1)


# Load annotations from Giladi's email  -----------------------------------


inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/tier3_annotation.txt"
assertthat::assert_that(file.exists(inf.meta))

dat.meta.orig <- fread(inf.meta)

jsplit <- lapply(split(dat.meta.orig, dat.meta.orig$marker), function(jdat) nrow(jdat))

dat.meta <- dat.meta.orig


# filter markers twice, one to remove redundant celltypes, second to keep only the top four named in a meaningful way 
if (!WithVpreb){
  markers.keep <- c("Car1", "core", "Siglech", "Prg2", "Ccl5", "Prss34", "Cd74", "Fcrla", "Ltf")
} else {
  markers.keep <- c("Car1", "core", "Siglech", "Prg2", "Ccl5", "Prss34", "Cd74", "Fcrla", "Ltf", "Vpreb1")  # add Vpreb1?
}
dat.meta <- subset(dat.meta, marker %in% markers.keep)

m2c <- MarkerToCelltype()
dat.meta$marker <- sapply(dat.meta$marker, function(x) ifelse(!is.null(m2c[[x]]), m2c[[x]], x))
print(unique(dat.meta$marker))

markers.keep.again <- c("Erythroblast", "HSCs", "Bcell", "Neutrophil")

dat.meta <- subset(dat.meta, marker %in% markers.keep.again)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.meta, aes(x = x, y = y, color = marker)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load data ---------------------------------------------------------------

indir <- "~/data/from_cluster/public_data/Giladi_et_al_2018"
assertthat::assert_that(dir.exists(indir))


inf.meta <- file.path(indir, "GSE92575_metadata.txt")

infs <- list.files(path = indir, pattern = "*.txt.gz", full.names = TRUE)
names(infs) <- sapply(infs, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

# get gene list
genes <- ReadGiladi(infs[[1]], remove.first.col = FALSE)$gene

dats <- lapply(infs, ReadGiladi, remove.first.col = TRUE) %>%
  bind_cols()

meta <- as.data.frame(fread(inf.meta, skip = 14))
rownames(meta) <- meta$well


# cells <- meta$well
cells <- dat.meta$V1

# integrate tier3 annotation into thiis
dats.filt <- dats[, ..cells]

dats.filt <- as.matrix(dats.filt)
rownames(dats.filt) <- genes

# split into tiers and then calculate entropy?


# maxcells for down sampling 


cells.by.marker <- lapply(split(x = dat.meta, f = dat.meta$marker), function(jsub){
  return(jsub$V1)
})

lapply(cells.by.marker, length)

ncells.by.marker <- lapply(cells.by.marker, length)

print(ncells.by.marker)

(maxcells <- min(unlist(ncells.by.marker)))





# keep only Bcells, eryth, HSPCs, and granulocytes



pbmc <- CreateSeuratObject(counts = dats.filt, project = "pbmc3k", min.cells = 3, min.features = 200, meta.data = meta)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))


all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", group.by = "tier") + scale_color_viridis_d()
DimPlot(pbmc, reduction = "umap", group.by = "tier", split.by = "tier", ncol = 5) + scale_color_viridis_d()


# Plot gene counts --------------------------------------------------------

stemcell.genes <- as.character(subset(markers, gene.module == "Stem genes")$gene)

jmodules <- c("Stem genes", "Monocyte", "Stage I neutrophil", "Stage II neutrophil")
marker.genes <- as.character(subset(markers, gene.module %in% jmodules)$gene)
marker.genes <- as.character(subset(markers, gene.module %in% jmodules)$gene)

# jgenes <- "Hlf"
# jgenes <- stemcell.genes[1:10]

jgenes <- "S100a8"
jgenes <- "Hlf"

# RidgePlot(pbmc, stemcell.genes, group.by = "tier")
# VlnPlot(pbmc, stemcell.genes, group.by = "tier")

FeaturePlot(pbmc, features = jgenes, order = TRUE)

# plot entropy in single cells


S.vec <- apply(pbmc@assays$RNA@counts, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
jmeta <- data.frame(S = S.vec)
rownames(jmeta) <- names(S.vec)
pbmc <- AddMetaData(object = pbmc, metadata = jmeta, col.name = "entropy")

FeaturePlot(pbmc, features = 'entropy') + scale_color_viridis_c(direction = 1)

# add real metadata???
dat.meta.celltypes <- as.data.frame(dat.meta)
rownames(dat.meta.celltypes) <- dat.meta.celltypes$V1

pbmc <- AddMetaData(object = pbmc, metadata = dat.meta.celltypes)

Seurat::DimPlot(pbmc, group.by = "marker")

if (make.plots){
  dev.off()
}

# Get markers -------------------------------------------------------------

Idents(pbmc) <- pbmc$marker
# de.output <- FindMarkers(pbmc, ident.1 = "Car1")

jmarkers <- unique(dat.meta.celltypes$marker)
names(jmarkers) <- jmarkers


de.output.lst <- FindAllMarkers(pbmc, only.pos = FALSE, logfc.threshold = 0, min.pct = 0, max.cells.per.ident = maxcells, random.seed = 0, verbose = TRUE, test.use = jmeth)

saveRDS(de.output.lst, file = outrds)

print(Sys.time() - jstart)
