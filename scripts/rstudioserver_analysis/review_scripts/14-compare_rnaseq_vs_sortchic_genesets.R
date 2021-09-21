# Jake Yeung
# Date of Creation: 2021-09-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/14-compare_rnaseq_vs_sortchic_genesets.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load scRNA-seq ----------------------------------------------------------

inf.baccin <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Baccin_2019/41556_2019_439_MOESM4_ESM/NicheData10x.rda")
load(inf.baccin, v=T)

print(NicheData10x)

# check UMAP
NicheData10x <- RunUMAP(NicheData10x, dims = 1:20)

# jcells.highlight <- 

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
 
DimPlot(NicheData10x, reduction = 'umap', group.by = 'metadata....experiment..') + 
  scale_color_manual(values = cbPalette) 

DimPlot(NicheData10x, reduction = 'umap')  + 
  scale_color_viridis_d()


# Remove mesenchymal cells  -----------------------------------------------

ctypes.mes <- c("Adipo-CAR", "Smooth muscle", "Osteo-CAR", "Osteoblasts", "Ng2+ MSCs", "Fibro/Chondro p.", "Stromal fibro.", "Arteriolar fibro.", "Endosteal fibro.", "Myofibroblasts", "Arteriolar ECs", "Sinusoidal ECs", "Schwann cells")
ctypes.mes.not <- NicheData10x@active.ident[which(!unique(NicheData10x@active.ident) %in% ctypes.mes)]

ctypes.immune <- c("Neutrophils", "pro-B", "small pre-B.", "large pre-B.", "B cell", "T cells", "NK cells", "Dendritic cells", "Monocytes")
ctypes.hspcs <- c("Neutro prog.", "Gran/Mono prog.", "Mono prog.", "LMPPs", "Eo/Baso prog.", "Erythroblasts", "Ery prog.", "Ery/Mk prog.")

ctypes.blood <- c(ctypes.immune, ctypes.hspcs)

# jcells.highlight <- WhichCells(object = NicheData10x, idents = ctypes.mes)
jcells.highlight <- WhichCells(object = NicheData10x, idents = c(ctypes.immune, ctypes.hspcs))
jcells.highlight.not <- WhichCells(object = NicheData10x, idents = ctypes.mes.not)

DimPlot(NicheData10x, reduction = 'umap', cells.highlight = jcells.highlight)  + 
  scale_color_viridis_d()

DimPlot(NicheData10x, reduction = 'umap', cells.highlight = jcells.highlight.not)  + 
  scale_color_viridis_d()


# Redo UMAP with only blood cells  ----------------------------------------

mat.blood <- NicheData10x@assays$RNA[, jcells.highlight]

Blood <- CreateSeuratObject(counts = mat.blood)
Blood <- NormalizeData(Blood)
Blood <- FindVariableFeatures(Blood)
Blood <- ScaleData(Blood)
Blood <- RunPCA(Blood, npcs = 30)
Blood <- RunUMAP(Blood, dims = 1:30)

ids <- Idents(NicheData10x)
jmeta <- ids[names(ids) %in% jcells.highlight]

Blood <- AddMetaData(Blood, metadata = jmeta, col.name = "celltype")

jout <- DimPlot(Blood, reduction = 'umap', group.by = "celltype") 

dat.umap <- jout$data %>%
  mutate(cell = rownames(jout$data))

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597", "#eb9d01", "#7fbedf")
m.umap <- ggplot(dat.umap, aes(x = UMAP_1, y = UMAP_2, color = celltype)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# # save output
# outobj <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq/Baccin_scRNAseq_bonemarrow_no_niche.", Sys.Date(), ".RData")
# save(Blood, dat.umap, file = outobj)
# 



# Load sortChIC data ------------------------------------------------------


inmain <- "jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq"
inf.tm.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_10kbTSS_LDA.rds")
inf.mat.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_10kbTSS_countmat.rds")
inf.meta.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_meta_withumap.txt")

tm.chic <- readRDS(inf.tm.chic)
mat.chic <- readRDS(inf.mat.chic)
dat.umap.chic.annot <- fread(inf.meta.chic)

m.chic <- ggplot(dat.umap.chic.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.imputed.chic <- t(log2(tm.chic$topics %*% tm.chic$terms))

# Color genesets onto UMAP  -----------------------------------------------

inf.gsets <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/post_submission/genesets/gene_sets_from_sortChIC.txt")
dat.gsets <- fread(inf.gsets)

jjset <- c("Eryths")

(jjsets <- unique(dat.gsets$jset))


outpdf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/post_submission/plots/scRNAseq_sortChIC_comparison.", Sys.Date(), ".pdf"))
pdf(outpdf, useDingbats = FALSE)
JFuncs::multiplot(m.umap, m.chic, cols = 2)
for(jjset in jjsets){
  
  # get mean exprs, plot onto UMAP 
  jgenes <- unique(subset(dat.gsets, jset == jjset)$symbol)
  jregions <- unique(subset(dat.gsets, jset == jjset)$gene)
  rows.keep <- rownames(mat.blood) %in% jgenes
  
  dat.imputed.chic.filt <- dat.imputed.chic[jregions, ]
  exprs.chic <- colMeans(dat.imputed.chic.filt)
  dat.exprs.chic <- data.frame(cell = names(exprs.chic), gexprs = exprs.chic, stringsAsFactors = FALSE) %>%
    left_join(., dat.umap.chic.annot)
  
  m.exprs.chic <- ggplot(dat.exprs.chic, aes(x = umap1, y = umap2, color = gexprs)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(paste0("+/-5kb TSS of set:", jjset), "sortChIC") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_viridis_c() 
  
  
  mat.sub <- mat.blood[rows.keep, ]
  dat.sub <- data.frame(cell = colnames(mat.sub), gexprs = colMeans(mat.sub), stringsAsFactors = FALSE) %>%
    left_join(., dat.umap)
  
  m.rnaseq <- ggplot(dat.sub, aes(x = UMAP_1, y = UMAP_2, color = gexprs)) + 
    geom_point() + 
    scale_color_viridis_c() + 
    theme_bw() + 
    ggtitle(paste0("Avg exprs of set:", jjset), "scRNAseq BM") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  # print(m.jset)
  JFuncs::multiplot(m.rnaseq, m.exprs.chic, cols = 2)
}
dev.off()


