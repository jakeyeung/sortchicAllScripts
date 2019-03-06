# Jake Yeung
# Date of Creation: 2019-02-22
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/downstream_MNN_analysis.R
# Downstream after running Seurat

library(Seurat)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

# Load data ---------------------------------------------------------------

bm.integrated <- readRDS("/tmp/bm_integrated.rds")
# bm.integrated <- readRDS("/tmp/bm_integrated.repressive.rds")

# Analyze -----------------------------------------------------------------

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = bm.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
bm.integrated <- ScaleData(object = bm.integrated, verbose = FALSE)
bm.integrated <- RunPCA(object = bm.integrated, npcs = 100, verbose = FALSE)
bm.integrated <- RunUMAP(object = bm.integrated, reduction = "pca", 
                         dims = 1:100)
p1 <- DimPlot(object = bm.integrated, reduction = "umap", group.by = "tech")
print(p1)

