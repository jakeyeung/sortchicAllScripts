# Jake Yeung
# Date of Creation: 2020-06-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/count_cells_per_cluster.R
# Count cells per cluster 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Indir -------------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/cell_to_cluster_tables_merged"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/cell_to_cluster_tables_merged.counts"
outf <- file.path("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/cell_to_cluster_tables_merged.counts", paste0("ncells_per_cluster.txt"))
dir.create(outdir)

cell.clusters <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(indir, paste0("BM_cell_to_clusters.", jmark, ".txt"))
  print(inf)
  clst <- fread(inf)
  return(clst)
})

cell.summary <- lapply(jmarks, function(jmark){
  jdat <- cell.clusters[[jmark]]
  jdat <- jdat %>% 
    group_by(cluster) %>%
    summarise(ncell = length(cell))
  jdat$mark <- jmark
  return(jdat)
})

cell.summary.long <- cell.summary %>% bind_rows()
  
fwrite(cell.summary.long, file = outf, quote = FALSE, sep = "\t")

