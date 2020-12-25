# Jake Yeung
# Date of Creation: 2020-12-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_heatmap_from_raw.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jmarks.act <- c("H3K4me1", "H3K4me3")

# Load data ---------------------------------------------------------------

inf.raw <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/raw_cuts_TSS_three_marks.rds"
dat.raw.lst <- readRDS(inf.raw)

# Load gene sets ----------------------------------------------------------

dat.genes <- GetGeneSets()

ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")

# Load metadata -----------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.", jmark, ".txt")
  dat.meta <- fread(inf)
  dat.meta$cluster <- factor(dat.meta$cluster, levels = ctypes)
  dat.meta <- dat.meta %>%
    arrange(cluster, jrep)
  return(dat.meta)
})


# Fix rownames ------------------------------------------------------------


rnames.act <- lapply(jmarks.act, function(jmark){
  rnames <- rownames(dat.raw.lst[[jmark]])
}) %>%
  unlist() %>%
  unique()  %>%
  gtools::mixedsort()

coords.act <- sapply(rnames.act, function(x) paste("chr", strsplit(x, split = ";")[[1]][[1]], sep = ""))

coord2rname <- hash::hash(coords.act, rnames.act)

rownames(dat.raw.lst$H3K27me3) <- sapply(rownames(dat.raw.lst$H3K27me3), function(x) AssignHash(x = x, jhash = coord2rname, null.fill = x))

# Get mat  ----------------------------------------------------------------

# K27me3



jmark <- "H3K27me3"
dat.raw.filt <- dat.raw.lst[[jmark]]
dat.raw.filt <- sweep(dat.raw.filt, MARGIN = 2, STATS = colSums(dat.raw.filt), FUN = "/")

genes.lst <- split(x = dat.genes$gene, f = dat.genes$jset)
clstrs.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)

pseudogene.mat <- SumAcrossClusters(t(dat.raw.filt), cnames.keep.lst = genes.lst)
pseudogene.mat <- do.call(rbind, pseudogene.mat)

# cells rearranged
cells.ordered <- dat.metas[[jmark]]$cell
pseudogenes.ordered <- ctypes



cols.hash.new <- hash::hash(dat.genes$jset, dat.genes$colorcode)

# Plot  -------------------------------------------------------------------



pseudogene.mat.reordered <- log2(pseudogene.mat[pseudogenes.ordered, cells.ordered] * 100 + 1)
pseudogene.mat.reordered <- t(apply(pseudogene.mat.reordered, 1, function(jrow) DescTools::Winsorize(jrow, probs = c(0.05, 0.99))))
pseudogene.mat.reordered <- apply(pseudogene.mat.reordered, 2, function(jcol) DescTools::Winsorize(jcol, probs = c(0.05, 0.99)))


# fix batch effects?
pseudogene.mat.reordered.long <- pseudogene.mat.reordered %>%
  data.table::melt()
colnames(pseudogene.mat.reordered.long) <- c("rname", "cell", "log2exprs")
plot(density(pseudogene.mat.reordered.long$log2exprs))


mat4hm.long.annot <- left_join(pseudogene.mat.reordered.long, dat.metas[[jmark]])
mat4hm.long.annot$jrep2 <- sapply(mat4hm.long.annot$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))

dat.adj <- mat4hm.long.annot %>%
  group_by(rname) %>%
  do(AdjustBatchEffect(.))

ggplot(dat.adj %>% filter(rname == "Granulocytes"), aes(x = cluster, y = log2exprsadj)) +
  geom_point() + 
  geom_boxplot() + 
  facet_wrap(~jrep)

mat.adj <- data.table::dcast(dat.adj, rname ~ cell, value.var = "log2exprsadj")
rownames(mat.adj) <- mat.adj$rname
mat.adj$rname <- NULL

mat.adj <- t(apply(mat.adj, 1, function(jrow) DescTools::Winsorize(jrow, probs = c(0.01, 0.99))))
mat.adj <- apply(mat.adj, 2, function(jcol) DescTools::Winsorize(jcol, probs = c(0.01, 0.99)))

heatmap3::heatmap3(pseudogene.mat.reordered, Rowv = NA, Colv = NA, ColSideColors = dat.metas[[jmark]]$colorcode, scale = "row", RowSideColors = sapply(rownames(pseudogene.mat), AssignHash, cols.hash.new), revC = TRUE, main = paste(jmark, "pseudogene"))

heatmap3::heatmap3(mat.adj[pseudogenes.ordered, cells.ordered], Rowv = NA, Colv = NA, ColSideColors = dat.metas[[jmark]]$colorcode, scale = "row", RowSideColors = sapply(rownames(pseudogene.mat), AssignHash, cols.hash.new), revC = TRUE, main = paste(jmark, "pseudogene"))

