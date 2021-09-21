# Jake Yeung
# Date of Creation: 2021-09-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/15-integrate_datasets_RNAseq_vs_sortChIC.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load RNAseq and sortChIC data  -------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inmain <- "jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq"

inf.rnablood <- file.path(hubprefix, inmain, "Baccin_scRNAseq_bonemarrow_no_niche.2021-09-05.RData")

load(inf.rnablood, v=T)  # Blood, dat.umap
dat.umap.rna.annot <- dat.umap %>%
  dplyr::rename(umap1 = UMAP_1,
                umap2 = UMAP_2)


inf.tm.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_10kbTSS_LDA.rds")
inf.mat.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_10kbTSS_countmat.rds")
inf.meta.chic <- file.path(hubprefix, inmain, "H3K4me1-sortChIC_Zeller_meta_withumap.txt")

tm.chic <- readRDS(inf.tm.chic)
mat.chic <- readRDS(inf.mat.chic)
dat.umap.chic.annot <- fread(inf.meta.chic)


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597",  "#ff9f7d", "#eb9d01", "#7fbedf")

m.chic <- ggplot(dat.umap.chic.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("sortChIC") + 
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.rna <- ggplot(dat.umap.rna.annot, aes(x = umap1, y = umap2, color = celltype)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("scRNAseq") + 
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.chic)
print(m.rna)


dat.imputed.chic <- t(log2(tm.chic$topics %*% tm.chic$terms))
dat.imputed.rna <- Blood@assays$RNA@counts


# Load gene sets ----------------------------------------------------------

inf.gsets <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/post_submission/genesets/gene_sets_from_sortChIC.txt")
dat.gsets <- fread(inf.gsets)
genes.keep <- unique(dat.gsets$symbol)
rna.genes.keep <- rownames(dat.imputed.rna)[rownames(dat.imputed.rna) %in% genes.keep]

set.seed(0)
dat.gsets.filt <- dat.gsets %>%
  filter(symbol %in% rna.genes.keep) %>%
  group_by(symbol) %>%
  # filter(rnk == min(rnk))
  slice(1)
  # top_n(x = ., n = -1, wt = rnk)

regions.keep.ordered <- dat.gsets.filt$gene
genes.keep.ordered <- dat.gsets.filt$symbol
print(length(regions.keep.ordered))

# Run CCA  ----------------------------------------------------------------

X <- as.matrix(dat.imputed.chic)[regions.keep.ordered, ]
Y <- as.matrix(dat.imputed.rna)[genes.keep.ordered, ]
rownames(X) <- genes.keep.ordered


assertthat::assert_that(nrow(X) == nrow(Y))

X <- scale(X, center = TRUE, scale = FALSE)
Y <- scale(Y, center = TRUE, scale = FALSE)

cca.results <- jCanonCor(X, Y, k = 20, l2.norm = TRUE)

cca.data <- rbind(cca.results$u, cca.results$v)

colnames(x = cca.data) <- paste0("CC", 1:20)
rownames(cca.data) <- c(colnames(X), colnames(Y))

# why does this matter? Maybe it does?
cca.data.flip <- apply(cca.data, MARGIN = 2, function(x){
  if(sign(x[1]) == -1) {
    x <- x * -1
  }
  return(x)
})


# how are these feature loadings calculated? Just projection
embeds <- cca.data.flip
# average across embeds
loads <- t(t(embeds) %*% rbind(t(X), t(Y)))
# take min function?
embeds1 <- cca.results$u
embeds2 <- cca.results$v

loads1 <- t(t(embeds1) %*% t(X))
loads2 <- t(t(embeds2) %*% t(Y))

loads.min <- matrix(mapply(SelectAbsMin, loads1, loads2), nrow = nrow(loads1), ncol = ncol(loads1), dimnames = list(rownames(loads1), colnames(loads1)))

cca1 <- 1
cca2 <- 2
plot(loads[, cca1], loads[, cca2], pch = 20)
# text(loads[, cca1], loads[, cca2],
#      labels = rownames(loads))
abline(v = 0, h = 0)


distfilt <- 60


loads.long <- data.frame(loads, bin = rownames(loads), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(dist = sqrt(CC1 ^ 2 + CC2 ^ 2)) %>%
  mutate(bin.lab = ifelse(dist > distfilt, bin, NA))



# Multipy loads -----------------------------------------------------------

X.proj <- t(loads) %*% X
Y.proj <- t(loads) %*% Y

XY <- cbind(X.proj, Y.proj)

colsvec <- c(rep("red", ncol(X.proj)), rep("blue", ncol(Y.proj)))
plot(XY[1, ], XY[2, ], col = colsvec)


# Seurat implementation ---------------------------------------------------


library(Seurat)

chic <- CreateSeuratObject(counts = X, meta.data = dat.umap.chic.annot, assay = "ChIC")
rna <- CreateSeuratObject(counts = Y, meta.data = dat.umap.rna.annot, assay = "RNA")

chic <- ScaleData(chic, do.scale = TRUE)
rna <- ScaleData(rna, do.scale = TRUE)

chic <- FindVariableFeatures(chic, nfeatures = 1000)
rna <- FindVariableFeatures(rna, nfeatures = 1000)

chic@meta.data$mark <- "H3K4me1"
rna@meta.data$mark <- "RNA"

marks.combined <- RunCCA(chic, rna, num.cc = 30)

p1 <- DimPlot(object = marks.combined, group.by = "mark", reduction = "cca", pt.size = 0.5, dims = 1:2)
p2 <- VlnPlot(object = marks.combined, features = "CC_1", group.by = "mark")



# Color by cell types  ------------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597", "#ff9f7d", "#eb9d01", "#7fbedf")


dat.cca <- marks.combined@reductions$cca@cell.embeddings
dat.umap.cca <- DoUmapAndLouvain(dat.cca, jsettings)



dat.umap.chic.annot.unif <- dat.umap.chic.annot %>% dplyr::select(c(cell, cluster)) %>%
  dplyr::mutate(experi = "ChIC")

dat.umap.rna.annot.unif <- dat.umap.rna.annot %>% dplyr::select(c(cell, celltype)) %>%
  dplyr::mutate(experi = "RNA") %>%
  dplyr::rename(cluster = celltype)

dat.umap.annot.unif <- rbind(dat.umap.chic.annot.unif, dat.umap.rna.annot.unif)


# color by experiment

dat.umap.cca.rna.chic.annot <- dat.umap.cca %>%
  left_join(., dat.umap.annot.unif, by = "cell") 

m.all <- ggplot(dat.umap.cca.rna.chic.annot, aes(x = umap1, y = umap2, color = experi)) + 
  geom_point(alpha = 0.5) + 
  ggtitle("Colored by dataset") + 
  xlab("CCA-UMAP1") + ylab("CCA-UMAP2") + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.chic <- ggplot(dat.umap.cca.rna.chic.annot %>% filter(experi == "ChIC"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(alpha = 1) + 
  xlab("CCA-UMAP1") + ylab("CCA-UMAP2") + 
  ggtitle("Colored by ChIC") + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.rna <- ggplot(dat.umap.cca.rna.chic.annot %>% filter(experi == "RNA"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(alpha = 1) + 
  xlab("CCA-UMAP1") + ylab("CCA-UMAP2") + 
  ggtitle("Colored by rna") + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

pdf(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/plots/chic_rna_CCA_comparison.", Sys.Date(), ".pdf"), useDingbats = FALSE)
print(m.all)
print(m.chic)
print(m.rna)
multiplot(m.chic, m.rna, cols = 2)
dev.off()

# save objects
rna.chic.comparison <- marks.combined
outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/analysis_scrnaseq_atacseq/scRNAseqBaccin_sortChICZeller_CCA_comparison.RData"
save(rna.chic.comparison, dat.umap.cca.rna.chic.annot, file = outf)
