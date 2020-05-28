# Jake Yeung
# Date of Creation: 2020-04-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/get_DE_genes_zebrafish_Chloe.R
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

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/For_Jake"

jlogfc.threshold <- 0
jmin.pct <- 0
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_rstudio"
outdir <- paste0("/home/jyeung/data/from_rstudioserver/zebrafish.", Sys.Date())  # run over night
dir.create(outdir)
outrds <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.rds"))
outpdf <- file.path(outdir, paste0("diff_exprs_Chloe_seurat.full.pdf"))

pdf(outpdf, useDingbats = FALSE)

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


# Clculate entropy --------------------------------------------------------


# calcualte entropy
S.vec <- apply(mat, 2, function(jcell) CalculateEntropy(jcell, normalize.p = TRUE))
jmeta <- data.frame(S = S.vec)
rownames(jmeta) <- names(S.vec)

meta.merge <- left_join(meta, data.frame(rname = rownames(jmeta), jmeta)) %>%
  filter(!is.na(S))

m <- ggplot(meta.merge, aes(x = S, group = celltype)) + facet_wrap(~celltype, ncol = 1) + geom_density(fill = "lightblue") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Entropy")
print(m)



# Load obj ----------------------------------------------------------------



zf <- CreateSeuratObject(counts = mat, project = "zf", min.cells = 10, min.features = 400)


# plot1 <- FeatureScatter(zf, feature1 = "nCount_RNA", feature2 = "percent.mt")
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
jhash.entropy <- hash(rownames(jmeta), jmeta$S)
zf@meta.data$experi <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash[[x]], NA)
})
zf@meta.data$clusterid <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash2[[x]]), jhash2[[x]], NA)
})
zf@meta.data$celltype <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash[[x]]), jhash3[[x]], NA)
})
zf@meta.data$entropy <- sapply(zf@meta.data$cell, function(x){
  ifelse(!is.null(jhash.entropy[[x]]), jhash.entropy[[x]], NA)
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

Idents(zf) <- zf$celltype

print("Running find markers")
de.output.lst <- FindAllMarkers(zf, only.pos = FALSE, logfc.threshold = jlogfc.threshold, min.pct = jmin.pct)
print("done Running find markers")


saveRDS(de.output.lst, file = outrds)

print(Sys.time() - jstart)