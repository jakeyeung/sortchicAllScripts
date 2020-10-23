# Jake Yeung
# Date of Creation: 2020-10-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/2-check_fit_outputs.filter_cells.FromAllMergedUnionRows.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)


# Paths -------------------------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")

hubprefix <- "/home/jyeung/hub_oudenaarden"

projmain <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/projects_after_unmixing.H3K4me1xH3K9me3/SetupObjs_AllMerged_UnionRows"))

inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.Robj"
assertthat::assert_that(file.exists(inf.dbl))

# inf.input <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/mouse_spikein_BMround2all.dbl_common_rows_match_dbl/mouse_spikein_BMround2all.dbl_common_rows_match_dbl_clstr_by_louvain_H3K4me1xH3K9me3.removeNA_TRUE.RData"
inf.input <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_UnionRows/SetupObjs_AllMerged_UnionRows.clstr_by_louvain_H3K4me1xH3K9me3.removeNA_FALSE.RData"
# inf.output <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/mouse_spikein_BMround2all.dbl_common_rows_match_dbl/unmix_mouse_spikein_BMround2all.dbl_common_rows_match_dbl_clstr_by_louvain_H3K4me1xH3K9me3.removeNA_TRUE.RData"
inf.output <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/SetupObjs_AllMerged_UnionRows/unmix_SetupObjs_AllMerged_UnionRows.clstr_by_louvain_H3K4me1xH3K9me3.removeNA_FALSE.RData"


assertthat::assert_that(file.exists(inf.input))
assertthat::assert_that(file.exists(inf.output))

indir.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_UnionRows/SetupObjs_AllMerged_UnionRows.clstr_by_louvain_H3K4me1xH3K9me3.removeNA_FALSE.RData"
# indir.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl"
assertthat::assert_that(dir.exists(indir.annot))
# inf.annot.k9me3 <- file.path(indir.annot, "cluster_tables_H3K9me3_BM_all_round2.txt")
# inf.annot.k4me1 <- file.path(indir.annot, "cluster_tables_H3K4me1_BM_all_round2.txt")

load(indir.annot, v=T)


dat.annot.k4me1 <- dat.louv$H3K4me1
dat.annot.k9me3 <- dat.louv$H3K9me3


# Load LDAs projections ---------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("H3K4me1", "H3K9me3")
names(jmarks) <- jmarks
# jmark <- "H3K4me1"

dat.umaps.lst <- lapply(jmarks, function(jmark){
  
  inf.projs <- file.path(projmain, paste0("project_unmixed_", jmark, ".RData"))
  load(inf.projs, v=T)
  
  tm.result <- posterior(out.objs$out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings) %>%
    rowwise() %>%
    mutate(stain = "single")
  
  tm.result.dbl <- AddTopicToTmResult(out.lda.predict)
  
  dat.umap.dbl <- DoUmapAndLouvain(tm.result.dbl$topics, jsettings = jsettings) %>%
    rowwise() %>%
    mutate(stain = "dbl")
  
  
  # Merge on one umap -------------------------------------------------------
  
  umap.out.orig <- umap(tm.result$topics, config = jsettings)
  umap.out.pred <- predict(umap.out.orig, data = tm.result.dbl$topics)
  
  dat.umap.out.orig <- data.frame(cell = rownames(umap.out.orig$layout), umap1 = umap.out.orig$layout[, 1], umap2 = umap.out.orig$layout[, 2], stain = "single", stringsAsFactors = FALSE)
  dat.umap.out.pred <- data.frame(cell = rownames(umap.out.pred), umap1 = umap.out.pred[, 1], umap2 = umap.out.pred[, 2], stain = "dbl", stringsAsFactors = FALSE)
  
  dat.umap.out.merge <- bind_rows(dat.umap.out.orig, dat.umap.out.pred)
  dat.umap.out.merge$mark <- jmark
  return(dat.umap.out.merge)
})


for (jmark in jmarks){
  m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = stain)) + 
    ggtitle(jmark) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~stain)
  print(m)
}



# Load UMAP from double stain ---------------------------------------------


load(inf.dbl, v=T)

tm.result.dbl <- posterior(out.lda)
tm.result.dbl <- AddTopicToTmResult(tm.result.dbl)

dat.umap.dbl <- DoUmapAndLouvain(topics.mat = tm.result.dbl$topics, jsettings = jsettings)


# Get assignments by clusters ---------------------------------------------


# inf.input <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/double_staining_input/", jprefix.input, "/", jprefix, "_", jmarks.dbl, ".removeNA_TRUE.RData")
# inf.output <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/double_staining_output/", jprefix.output, "/unmix_", jprefix, "_", jmarks.dbl, ".removeNA_TRUE.RData")


load(inf.input, v=T)
load(inf.output, v=T)



# Make grid ---------------------------------------------------------------


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


# coords.dbl.annots <- left_join(coords.dbl, annots.dat)
coords.dbl.annots <- coords.dbl

dat.umap.dbl.merge <- left_join(dat.umap.dbl, coords.dbl.annots) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

m.check.repress <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85")  + facet_wrap(~plate)

m.check.act <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = louv.act)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey85")  + facet_wrap(~plate)

m.check.w <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = w)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + 
  facet_wrap(~plate)

print(m.check.repress)
print(m.check.act)
print(m.check.w)

# coords.dbl.annots$louv.act <- factor(coords.dbl.annots$louv.act, levels = louv.act.ordering)
# coords.dbl.annots$louv.repress <- factor(coords.dbl.annots$louv.repress, levels = louv.repress.ordering)


library(ggforce)
m.grid <- ggplot(coords.dbl.annots, aes(x = louv.act, y = louv.repress, louv.act)) + 
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) + 
  theme_bw() + 
  theme(aspect.ratio=0.6) + 
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) + 
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned") 

print(m.grid)

m.grid.w <- ggplot(coords.dbl.annots, aes(x = louv.act, y = louv.repress, col = w)) + 
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) + 
  xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) + 
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned") 
print(m.grid.w)

# Rename  -----------------------------------------------------------------


dat.annot.k9me3 <- fread(inf.annot.k9me3) %>%
  mutate(louvain = paste("louvain", louvain, sep = ""))

m.k9me3 <- ggplot(dat.annot.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

dat.annot.k4me1 <- fread(inf.annot.k4me1) %>%
  mutate(louvain = paste("louvain", louvain, sep = ""))

m.k4me1 <- ggplot(dat.annot.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


JFuncs::multiplot(m.grid, m.k9me3, cols = 2)




# Plot UMAP H3K4me1 single and double -------------------------------------

jmark <- "H3K4me1"
m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = stain)) + 
  ggtitle(jmark) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~stain)
print(m)

# label from dat.annot
dat.umap.k4me1 <- left_join(dat.umaps.lst[["H3K4me1"]], subset(dat.annot.k4me1, select = c(cell, cluster))) %>%
  left_join(., coords.dbl.annots) %>%
  rowwise() %>%
  mutate(cluster = ifelse(is.na(cluster), louv.act, cluster))

dat.umap.k9me3 <- left_join(dat.umaps.lst[["H3K9me3"]], subset(dat.annot.k9me3, select = c(cell, cluster))) %>%
  left_join(., coords.dbl.annots) %>%
  rowwise() 


m <- ggplot(dat.umap.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  ggtitle("H3K4me1") + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~stain)
print(m)

m <- ggplot(dat.umap.k4me1, aes(x = umap1, y = umap2, color = louv.act)) + 
  ggtitle(jmark) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~stain)
print(m)

m <- ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  ggtitle("H3K9me3") + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~stain)
print(m)

m <- ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = louv.act)) + 
  ggtitle("H3K9me3") + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~stain)
print(m)


# Check if there are double cells that are bad ----------------------------


print(m.check.act)
print(m.check.repress)


inf.annot.dbl <- file.path(indir.annot, "cluster_tables_H3K4me1xH3K9me3_BM_all_round2.txt")
dat.annot.dbl <- fread(inf.annot.dbl)

m.grid <- ggplot(coords.dbl.annots, aes(x = louv.act, y = louv.repress, louv.act)) + 
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) + 
  theme_bw() + 
  theme(aspect.ratio=0.6) + 
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) + 
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned") 

print(m.grid)





# 
# 
# m.grid <- ggplot(coords.dbl.annots %>% rowwise() %>% mutate(louv.act = AssignHash(louv.act, clusterhash, null.fill = NA)), aes(x = louv.act, y = louv.repress, louv.act)) + 
#   geom_point(position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) + 
#   theme_bw() + 
#   theme(aspect.ratio=0.35) + 
#   scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) + 
#   ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned") 
# 
# print(m.grid)
