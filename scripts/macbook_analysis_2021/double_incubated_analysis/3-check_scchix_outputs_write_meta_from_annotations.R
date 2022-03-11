# Jake Yeung
# Date of Creation: 2022-01-04
# File: ~/projects/scchic/scripts/macbook_analysis_2021/double_incubated_analysis/3-check_scchix_outputs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/from_cluster_2021/snakemake_H3K4me1_H3K9me3_again/snakemake_H3K4me1_H3K9me3_wlimits/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_H3K4me1-H3K9me3.RData"
load(inf, v=T)

fits.out <- act.repress.coord.lst

w.lst <- lapply(fits.out, function(jout)jout$w)

plot(w.lst %>% unlist())



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


m.grid.w.filt <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w, lnprob == 0)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_viridis_c() + 
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.w.filt)


# Load metadata  ----------------------------------------------------------

inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/from_cluster_2021/double_stain_outputs/scchix_outputs.H3K4me1-H3K9me3.2021-02-11.setseed.RData"
load(inf.meta, v=T)



coords.dbl.annot <- left_join(coords.dbl, subset(dat.umap.long.annot, select = c(cell, louv.repress.impute, louv.act.impute)))

m.grid <- ggplot(coords.dbl.annot, aes(x = louv.act, y = louv.repress, color = louv.act.impute)) +
  geom_point(alpha = 1, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)


jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks
jmarks.all <- c("H3K4me1", "H3K9me3", "H3K4me1-H3K9me3"); names(jmarks.all) <- jmarks.all
inf.meta.singles.lst <- lapply(jmarks, function(jmark){
  inf.meta.singles <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/metadata/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  return(inf.meta.singles)
})

dat.meta.singles.lst <- lapply(inf.meta.singles.lst, function(jinf){
  fread(jinf)
})

dat.meta.doubles <- dat.umap.long.annot %>%
  dplyr::select(cell, umap1, umap2, batch, louv.act.impute) %>%
  dplyr::rename(cluster = louv.act.impute)

dat.meta.doubles.lst <- list("H3K4me1-H3K9me3" = dat.meta.doubles)

# Make prety metadatas  ----------------------------------------------------

# write cell, umap1, umap2, and cluster

dat.metas.all.lst <- c(dat.meta.singles.lst, dat.meta.doubles.lst)

dat.metas.all.clean.lst <- lapply(dat.metas.all.lst, function(x){
  subset(x, select = c(cell, umap1, umap2, batch, cluster))
})



# Check mats --------------------------------------------------------------

# indir.mats <- "/Users/yeung/data/scchic/from_cluster_2021/metadata/metadata_for_scchix/snakemake_H3K4me1_H3K9me3.objs_from_LDA"
# infs.mats <- lapply(jmarks.all, function(jmark){
#   file.path(indir.mats, paste0("countmat_output_filt.", jmark, ".rds"))
# })
# 
# count.mats.lst <- lapply(infs.mats, function(jinf) readRDS(jinf))

# count.mats.filt.lst <- lapply(jmarks.all, function(jmark){
#   jmat <- count.mats.lst[[jmark]]
#   cells.keep <- dat.metas.all.clean.lst[[jmark]]$cell
#   cells.keep.i <- colnames(jmat) %in% cells.keep
#   jmat.filt <- jmat[, cells.keep.i]
#   jmat.filt <- jmat.filt[, !duplicated(colnames(jmat.filt))]
#   assertthat::assert_that(ncol(jmat.filt) > 0)
#   # remove duplicates
#   return(jmat.filt)
# })


# save ooutput
outdir <- "/Users/yeung/data/scchic/from_cluster_2021/metadata/metadata_for_scchix/snakemake_H3K4me1_H3K9me3.objs_from_LDA.cluster_from_annotations.ready_to_copy"
for (jmark in jmarks.all){
  print(jmark)
  # jmat <- count.mats.lst[[jmark]]
  jmeta <- dat.metas.all.clean.lst[[jmark]]
  # saveRDS(object = jmat, file.path(outdir, paste0("countmat_output_filt.", jmark, ".rds")))
  saveRDS(object = jmeta, file.path(outdir, paste0("celltyping_output_filt.", jmark, ".rds")))
}


## # check TM
## load("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/metadata/metadata_for_scchix/snakemake_H3K4me1_H3K9me3.objs_from_LDA/lda_output_filt.H3K4me1.rds", v=T)
