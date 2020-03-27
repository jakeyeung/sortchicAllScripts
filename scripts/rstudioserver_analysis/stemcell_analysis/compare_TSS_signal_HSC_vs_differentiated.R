# Jake Yeung
# Date of Creation: 2020-03-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/compare_TSS_signal_HSC_vs_differentiated.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Functions ---------------------------------------------------------------


MakeLongAndMerge <- function(mat, dat.merge, count.cname = "count.mark1", shorten.genenames = TRUE){
  # mat is sparse matrix genes in rows, cells in columns.
  # dat.merge are annotations with colnames in "cell" matching rownames in mat
  if (shorten.genenames){
    rownames(mat) <- sapply(rownames(mat), function(x) strsplit(x, ";")[[1]][[2]])
  }
  dat <- data.table::melt(as.matrix(mat))
  colnames(dat) <- c("gene", "cell", count.cname)
  # dat$gene <- sapply(dat$gene.full, function(x) strsplit(x, ";")[[1]][[2]])
  dat.merge <- left_join(dat, dat.merge)
  return(dat.merge)
}

# Load data ---------------------------------------------------------------


jmark.act <- "H3K4me3"
jmark.repress <- "H3K27me3"

inf.act <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/proxfrac.TSS_10000.", jmark.act, ".RData")
inf.repress <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/proxfrac.TSS_10000.", jmark.repress, ".RData")

# dat.de is common across marks 

load(inf.act, v=T)
mat.act <- mat
dat.merge.act <- dat.merge

load(inf.repress, v=T)
mat.repress <- mat
dat.merge.repress <- dat.merge

# Load DE genes -----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
dat.de <- readRDS(inf.de)

unique(dat.de$cluster)


# Make long dat with active and repress in one row  -----------------------

dat.act <- MakeLongAndMerge(mat.act, dat.merge.act, count.cname = "counts.act", shorten.genenames = TRUE)
dat.repress <- MakeLongAndMerge(mat.repress, dat.merge.repress, count.cname = "count.repress", shorten.genenames = TRUE)

# Sanity check: houskeeping and neutro genes ------------------------------

print(unique(dat.de$cluster))

jclst <- "Fcrla"
jclst <- "Ccl5"
jclst <- "Siglech"
jgenes <- subset(dat.de, cluster == jclst & p_val_adj < 0.01 & avg_logFC > 0)$gene

# neutro genes
jsub <- subset(dat.act, gene %in% jgenes) %>%
  group_by(cell) %>%
  summarise(counts.act = sum(counts.act)) %>%
  left_join(., dat.merge.act) %>%
  mutate(count.norm = counts.act / total.cuts)

ggplot(jsub, aes(x = umap1, y = umap2, color = log2(count.norm * 10^3 + 1))) + 
  scale_color_viridis_c() + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jclst)


# To do zscore, we need to first pseudobulk  ------------------------------

dat.act.bulk <- dat.act %>%
  group_by(cluster, gene) %>%
  summarise(counts.act = sum(counts.act))

# To merge ACT and REPRESS, we need to annotate cells by cluster  ---------



