# Jake Yeung
# Date of Creation: 2019-04-01
# File: ~/projects/scchic/scripts/scripts_analysis/make_tables/write_umap_tables_for_spring.R
# Write tables for spring exploration


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

for (jmark in jmarks){
  cell.topic.mat <- tm.result.lst[[jmark]]$topics
  
  imputed.mat <- as.data.frame(tm.result.lst[[jmark]]$topics %*% tm.result.lst[[jmark]]$terms)
  
  cell.bcs <- rownames(imputed.mat)
  
  imputed.mat.df <- data.frame(barcodes = cell.bcs, imputed.mat)
  
  print(dim(imputed.mat.df))
  
  # write imputed.mat to output 
  outf <- paste0("/Users/yeung/data/scchic/tables/imputed_matrix/", jmark, "_imputed_matrix_transpose.txt")
  outf.genelst <- paste0("/Users/yeung/data/scchic/tables/imputed_matrix/", jmark, "_genelist.txt")
  
  data.table::fwrite(imputed.mat.df, file = outf)
  system(paste0("gzip ", outf))
  
  # create gene list
  gene.lst <- colnames(imputed.mat)
  write(gene.lst, file = outf.genelst)
}

