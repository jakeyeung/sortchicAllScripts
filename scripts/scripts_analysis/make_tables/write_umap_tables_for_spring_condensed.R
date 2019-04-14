# Jake Yeung
# Date of Creation: 2019-04-08
# File: ~/projects/scchic/scripts/scripts_analysis/make_tables/write_umap_tables_for_spring_condensed.R
# Write just the UMAP matrix

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.RData"

load(inf, v=T)


# Write topic models list  ------------------------------------------------

# do only one mark
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
# jmark <- "H3K4me1"

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")

for (jmark in jmarks){
  cell.topic.mat <- tm.result.lst[[jmark]]$topics
  cell.bcs <- rownames(cell.topic.mat)
  colnames(cell.topic.mat) <- paste("Topic", colnames(cell.topic.mat), sep = "_")
  imputed.mat.df <- data.frame(barcodes = cell.bcs, cell.topic.mat)
  
  # write imputed.mat to output 
  outf <- paste0("/Users/yeung/data/scchic/tables/imputed_matrix/", jmark, "_imputed_matrix_transpose_cell_loadings.txt")
  outf.genelst <- paste0("/Users/yeung/data/scchic/tables/imputed_matrix/", jmark, "_genelist_cell_loadings.txt")
  
  data.table::fwrite(imputed.mat.df, file = outf)
  system(paste0("gzip ", outf))
  
  # create gene list
  gene.lst <- colnames(cell.topic.mat)
  write(gene.lst, file = outf.genelst)
  
  # create cell groupings 
  dat.groupings <- subset(dat.umap.long.new.lst[[jmark]], select = c(cell, louvain))
  dat.groupings$colorlabel <- cbPalette[as.numeric(dat.groupings$louvain)]
  data.table::fwrite(dat.groupings, file = paste0("/Users/yeung/data/scchic/tables/imputed_matrix/", jmark, "_color_labels.txt"))
}

