# Jake Yeung
# Date of Creation: 2020-06-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/check_exprs_of_genes_across_pseudobulk_rnaseq.R
# Load RNAseq 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds"

pdf(outpdf, file = paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds/exprs_mixed_states_BM_WKM.", Sys.Date(), ".pdf"), useDingbats = FALSE)


# Load bins  --------------------------------------------------------------

inf.mixed.bm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds/fit_poisson_model_on_TSS.Downstream2DClouds.2020-06-14.WithCI.MixedBins.txt"
inf.mixed.wkm <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds/fit_poisson_model_on_TSS.Downstream2DClouds.ZF.2020-06-14.ClusterRenamed.WithCI.MixedBin.txt"

dat.mixed.bins.bm <- fread(inf.mixed.bm)
dat.mixed.bins.wkm <- fread(inf.mixed.wkm)

mixed.ens.bm <- unique(dat.mixed.bins.bm$ens)
mixed.ens.wkm <- unique(dat.mixed.bins.wkm$ens)

# Load RNAseq  ------------------------------------------------------------

# mouse
inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

# zebrafish
inf.zf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.get_DE_genes_from_pbulk_scrnaseq/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.exprsmax_5.logfcmin_2.logfcmax_1.RData"
load(inf.zf.de, v=T)


# Plot BM mixed sttate gene expression across states ----------------------

ggplot(dat.sorted.norm.long, aes(x = exprs, fill = CellType)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")

ggplot(dat.sorted.norm.long, aes(x = exprs)) + geom_density(alpha = 0.25, fill = 'blue') + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")

ggplot(dat.sorted.norm.long %>% filter(gene %in% mixed.ens.bm), aes(x = CellType, y = exprs)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")

ggplot(dat.sorted.norm.long %>% filter(gene %in% mixed.ens.bm), aes(x = CellType, y = logFC)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")

ggplot(dat.sorted.norm.long %>% filter(gene %in% mixed.ens.bm), aes(x = CellType, y = zscore)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")


ggplot(dat.giladi.norm.long, aes(x = exprs, fill = pbulk)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Mouse BM Giladi scRNAseq pseudobulk")
  

ggplot(dat.giladi.norm.long, aes(x = exprs)) + geom_density(alpha = 0.25, fill = 'blue') + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Mouse BM Giladi scRNAseq pseudobulk")

ggplot(dat.giladi.norm.long %>% filter(ens %in% mixed.ens.bm), aes(x = pbulk, y = exprs)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM Giladi scRNAseq pseudobulk")

ggplot(dat.giladi.norm.long %>% filter(ens %in% mixed.ens.bm), aes(x = pbulk, y = logfc)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM Giladi scRNAseq pseudobulk")

ggplot(dat.giladi.norm.long %>% filter(ens %in% mixed.ens.bm), aes(x = pbulk, y = zscore)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM Giladi scRNAseq pseudobulk")


# Plot WKM mixed sttate gene expression across states ----------------------


ggplot(pbulk.zf.ctypefilt.long, aes(x = log2p1counts, fill = pbulk)) + geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("ZF WKM Baron scRNAseq pseudobulk")

ggplot(pbulk.zf.ctypefilt.long, aes(x = log2p1counts)) + geom_density(alpha = 0.25, fill = 'blue') + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("ZF WKM Baron scRNAseq pseudobulk")

ggplot(pbulk.zf.ctypefilt.long %>% filter(ens %in% mixed.ens.wkm), aes(x = pbulk, y = log2p1counts, fill = pbulk)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("ZF WKM Baron scRNAseq pseudobulk")

ggplot(pbulk.zf.ctypefilt.long %>% filter(ens %in% mixed.ens.wkm), aes(x = pbulk, y = log2fc, fill = pbulk)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("ZF WKM Baron scRNAseq pseudobulk")

ggplot(pbulk.zf.ctypefilt.long %>% filter(ens %in% mixed.ens.wkm), aes(x = pbulk, y = zscore, fill = pbulk)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("ZF WKM Baron scRNAseq pseudobulk")

dev.off()


# Overlap HSCs with mixed states ------------------------------------------

hsc.mixed.overlap <- intersect(de.ens.sorted.stringent$HSCs, mixed.ens.bm)

jens.vec <- unique(dat.sorted.norm.long$gene)
jgenes.vec <- JFuncs::EnsemblGene2Gene(jens.vec, return.original=TRUE)
g2e.mouse <- hash::hash(as.character(jgenes.vec), as.character(jens.vec))

# plot some genes for fun 
jgene <- "F3"
genes.keep <- c("F3", "Rflnb", "Olfm4", "Prom1", "Nucb2", "Ddx60", "Celsr3", "Fcer1a", "Mogat2", "Sgms2", "Edrna", "Pilrb2", "Aldh3b2", "S100a7a", "Chst13", "Lgals4", "Tmem45a2", "Pcp4", "Rnd1", "Map3k15", "Atrn", "Rps6ka2", "Prok2", "Kcnn7", "Ifnlr1")
print(length(unique(genes.keep)))
# jens <- g2e.hash[[jgene]]
ens.keep <- sapply(genes.keep, AssignHash, g2e.mouse)

ggplot(dat.sorted.norm.long %>% filter(gene %in% ens.keep), aes(x = CellType, y = zscore)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")

ggplot(dat.sorted.norm.long %>% filter(gene %in% ens.keep), aes(x = CellType, y = logFC)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")

ggplot(dat.sorted.norm.long %>% filter(gene %in% ens.keep), aes(x = CellType, y = exprs)) + 
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Mouse BM sorted celltypes bulk RNAseq")


