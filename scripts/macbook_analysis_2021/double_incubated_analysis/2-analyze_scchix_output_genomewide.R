# Jake Yeung
# Date of Creation: 2022-01-03
# File: ~/projects/scchic/scripts/macbook_analysis_2021/double_incubated_analysis/2-analyze_scchix_output_genomewide.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

inf <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.50kb_genomewide/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.RData"
load(inf, v=T)

act.repress.coord.lst[[1]]


fits.out <- act.repress.coord.lst

w.lst <- sapply(fits.out, function(x) x$w)

# remove 0.01 or 0.99
cells.remove.i <- which(w.lst >= 0.99 | w.lst <= 0.01)
if (length(cells.remove.i) > 0){
  cells.remove <- names(w.lst)[cells.remove.i]
  fits.out[[cells.remove]] <- NULL
}


# if louvains are now from clusters need eto rethink jcoord
cell.vec <- names(fits.out)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- fits.out[[jcell]]
  jweight <- fits.out[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)
  
  # rows are active, columns are repress I THINK?
  # TODO: assumes underscores be careful!
  jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
  jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
  
  if (grepl("_", jlouv.act)){
    jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
  }
  if (grepl("_", jlouv.repress)){
    jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
  }
  
  out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()




library(ggforce)
jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = louv.act)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)


m.grid.w <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_viridis_c() + 
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.w)


m.grid.w.filt <- ggplot(coords.dbl %>% filter(lnprob == 0), aes(x = louv.act, y = louv.repress, color = w, lnprob == 0)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_viridis_c() + 
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.w.filt)


