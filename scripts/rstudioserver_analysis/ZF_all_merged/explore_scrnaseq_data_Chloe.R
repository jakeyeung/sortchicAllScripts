# Jake Yeung
# Date of Creation: 2020-04-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/explore_scrnaseq_data_Chloe.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(Seurat)

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset_PenalizedNegBinRegress.rds")

dat <- readRDS(inf)


# Load pseudobulk gene expression?  ---------------------------------------

inf.scrnaseqbulk <- "/home/jyeung/data/from_rstudioserver/zebrafish.poisson.SameNbrCells.2020-04-30/WKM_pseudobulk_scrnaseq_downsampled.2020-04-30.SameNbrCells.RData"
load(inf.scrnaseqbulk, v=T)

pbulk.ctypefilt.long <- pbulk.ctypefilt.long %>%
  rowwise() %>%
  mutate(ens = strsplit(as.character(gene), split = "_")[[1]][[1]])

jgenesfull <- rownames(mat.pbulk.ds.ctypefilt)
jens <- sapply(jgenesfull, function(x) strsplit(x, "_")[[1]][[1]], USE.NAMES = FALSE)
jgenes <- sapply(jgenesfull, function(x) strsplit(x, "_")[[1]][[2]], USE.NAMES = FALSE)

g2e <- hash(jgenes, jens)


# Check gene expression output --------------------------------------------


# load from Seurat
inf.seurat <- "/home/jyeung/data/from_rstudioserver/zebrafish.poisson.SameNbrCells.2020-04-30/diff_exprs_Chloe_seurat.full.ctypefilt.SameNbrCells.rds"
dat.seurat <- readRDS(inf.seurat) %>%
  rowwise() %>%
  mutate(ens = strsplit(gene, "-")[[1]][[1]])


jgenes.hspc <- c("meis1b", "pmp22b", "ahnak", "krt8", "anxa5b", "mrc1b", "pa2g4a", "npm1a", "ahcy", "adh5", "fabp3", "myb")
jens.keep <- sapply(jgenes.hspc, AssignHash, g2e, null.fill = NA)

# plot HSPC-specific genes
jsub.hspc <- subset(dat.seurat, cluster == "HSPCs" & p_val_adj < 0.001 & avg_logFC > 1) 

jens.keep <- jsub.hspc$ens

ggplot(pbulk.ctypefilt.long %>% filter(ens %in% jens.keep), aes(x = pbulk, y = log2fc)) + geom_boxplot()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pbulk.ctypefilt.long %>% filter(ens %in% jens.keep), aes(x = pbulk, y = log2p1counts)) + geom_boxplot()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(pbulk.ctypefilt.long %>% filter(ens %in% jens.keep), aes(x = pbulk, y = counts)) + geom_boxplot()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(pbulk.ctypefilt.long, aes(x = log2p1counts)) + geom_density()
