# Jake Yeung
# Date of Creation: 2020-10-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/4-match_UMAPs_K4me1_and_K9me3.KeepAllCellFewerRepressedClusters.R
# description

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

remove.na <- TRUE
joutsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters"
# joutsuffix <- "UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_scChIX_output/check_louvs_umaps.", joutsuffix, ".", Sys.Date(), ".AnnotateUMAPsByCelltype.pdf")
outfinal <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_scChIX_output/match_UMAP_assign_clusters_objs.", joutsuffix, ".", Sys.Date(), ".final.RData")

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#851663", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")

hubprefix <- "/home/jyeung/hub_oudenaarden"

projmain <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/projects_after_unmixing.H3K4me1xH3K9me3/SetupObjs_AllMerged_", joutsuffix))

inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.Robj"
assertthat::assert_that(file.exists(inf.dbl))

# inf.input <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/mouse_spikein_BMround2all.dbl_common_rows_match_dbl/mouse_spikein_BMround2all.dbl_common_rows_match_dbl_clstr_by_louvain_H3K4me1xH3K9me3.removeNA_TRUE.RData"
inf.input <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_", joutsuffix, "/SetupObjs_AllMerged_", joutsuffix, ".clstr_by_louvain_H3K4me1xH3K9me3.removeNA_", remove.na, ".RData")
# inf.output <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/mouse_spikein_BMround2all.dbl_common_rows_match_dbl/unmix_mouse_spikein_BMround2all.dbl_common_rows_match_dbl_clstr_by_louvain_H3K4me1xH3K9me3.removeNA_TRUE.RData"
inf.output <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/SetupObjs_AllMerged_", joutsuffix, "/unmix_SetupObjs_AllMerged_", joutsuffix, ".clstr_by_louvain_H3K4me1xH3K9me3.removeNA_", remove.na, ".RData")


assertthat::assert_that(file.exists(inf.input))
assertthat::assert_that(file.exists(inf.output))

inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/SetupObjs_AllMerged_", joutsuffix, "/SetupObjs_AllMerged_", joutsuffix, ".clstr_by_louvain_H3K4me1xH3K9me3.removeNA_", remove.na, ".RData")
# indir.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl"
# assertthat::assert_that(dir.exists(inf.annot))
assertthat::assert_that(file.exists(inf.annot))
# inf.annot.k9me3 <- file.path(indir.annot, "cluster_tables_H3K9me3_BM_all_round2.txt")
# inf.annot.k4me1 <- file.path(indir.annot, "cluster_tables_H3K4me1_BM_all_round2.txt")

load(inf.annot, v=T)


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

m.check.repress <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y = umap2, color = louv.repress)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
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


# dat.annot.k9me3 <- fread(inf.annot.k9me3) %>%
#   mutate(louvain = paste("louvain", louvain, sep = ""))

m.k9me3 <- ggplot(dat.annot.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# dat.annot.k4me1 <- fread(inf.annot.k4me1) %>%
#   mutate(louvain = paste("louvain", louvain, sep = ""))

m.k4me1 <- ggplot(dat.annot.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

print(m.k4me1)
print(m.k9me3)

JFuncs::multiplot(m.grid, m.k9me3, cols = 2)
JFuncs::multiplot(m.grid, m.k4me1, cols = 2)




# Plot UMAP H3K4me1 single and double -------------------------------------

jmark <- "H3K4me1"
m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = stain)) + 
  ggtitle(jmark) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
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
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~stain)
print(m)

m <- ggplot(dat.umap.k4me1, aes(x = umap1, y = umap2, color = louv.act)) + 
  ggtitle(jmark) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~stain)
print(m)

m <- ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  ggtitle("H3K9me3") + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~stain)
print(m)

m <- ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = louv.act)) + 
  ggtitle("H3K9me3") + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~stain)
print(m)


# Check if there are double cells that are bad ----------------------------


print(m.check.act)
print(m.check.repress)

JFuncs::multiplot(m.grid, m.k4me1, cols = 2)
JFuncs::multiplot(m.grid, m.k9me3, cols = 2)



# Connect the UMAPs  ------------------------------------------------------


# shift k4me1 left, shift k9me3 right
dat.umap.k4me1.shift <- dat.umap.k4me1 %>%
  ungroup() %>%
  mutate(umap1 = unlist(scale(umap1, center = TRUE, scale = TRUE)), 
         umap1.shift = umap1 - 2)

dat.umap.k9me3.shift <- dat.umap.k9me3 %>%
  ungroup() %>%
  mutate(umap1 = unlist(scale(umap1, center = TRUE, scale = TRUE)), 
         umap1.shift = umap1 + 2)

dat.umap.merged <- rbind(dat.umap.k4me1.shift, dat.umap.k9me3.shift)

# put everything in one UMAP 

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")

# check louvs active
jlouvs <- unique(dat.umap.merged$louv.act)
jlouvs.act <- jlouvs[!is.na(jlouvs)]

jlouvs <- unique(dat.umap.merged$louv.repress)
jlouvs.repress <- jlouvs[!is.na(jlouvs)]


pdf(file = outpdf, useDingbats = FALSE)
# do acts
m.k4me1.act <- ggplot(dat.annot.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle("H3K4me1 clusters") + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.k4me1.act)

m.k9me3.rep <- ggplot(dat.annot.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle("H3K9me3 clusters") + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.k9me3.rep)

for (jlouv in jlouvs.act){
  m <- ggplot(dat.umap.merged %>% filter(stain == "dbl") %>% mutate(louv.act = louv.act == jlouv), aes(x = umap1.shift, y = umap2, group = cell, color = louv.act)) + 
    geom_line(alpha = 0.01) + 
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.25, size = 1) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(paste("Active:", jlouv)) + 
    theme_bw() + theme(aspect.ratio=0.75, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


m.k9me3.rep <- ggplot(dat.annot.k9me3, aes(x = umap1, y = umap2, color = cluster)) + 
  scale_color_manual(values = cbPalette) + 
  ggtitle("H3K9me3 Repress clusters") + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.k9me3.rep)

for (jlouv in jlouvs.repress){
  m <- ggplot(dat.umap.merged %>% filter(stain == "dbl") %>% mutate(louv.repress = louv.repress == jlouv), aes(x = umap1.shift, y = umap2, group = cell, color = louv.repress)) + 
    geom_line(alpha = 0.01) + 
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.25, size = 1) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(paste("Repres:", jlouv)) + 
    theme_bw() + theme(aspect.ratio=0.75, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


# Check var ---------------------------------------------------------------


inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.Robj"
load(inf.dbl, v=T)
tm.result.dbl <- posterior(out.lda)

dat.impute.log <- t(log2(tm.result.dbl$topics %*% tm.result.dbl$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.merged.var <- left_join(dat.umap.merged, dat.var)

ggplot(dat.umap.merged.var %>% filter(mark == "H3K4me1"), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)


ggplot(dat.umap.merged.var %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# filter out low var?
varfilt <- 1
ggplot(dat.umap.merged.var %>% filter(mark == "H3K9me3" & cell.var.within.sum.norm >= varfilt), aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.umap.merged.var %>% filter(mark == "H3K4me1" & cell.var.within.sum.norm >= varfilt), aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.umap.merged.var %>% filter(mark == "H3K9me3" & cell.var.within.sum.norm >= varfilt), aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.umap.merged.var %>% filter(mark == "H3K9me3" & (cell.var.within.sum.norm >= varfilt | umap2 < -2)), 
       aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

ggplot(dat.umap.merged.var %>% filter(mark == "H3K9me3" & cell.var.within.sum.norm >= 0) %>% mutate(louv.act = ifelse(louv.act == "louvain5x1", "IsStem", "NotStem")), aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)


# Annotate louv.act  ------------------------------------------------------
# 

# 
# merge.louvs.active.hash <- list("louvain15" = "Eryths",
#                                 "louvain2" = "Eryths",
#                                 "louvain13" = "DCs",
#                                 "louvain6" = "DCs",
#                                 "louvain7" = "Granus",
#                                 "louvain4" = "Granus",
#                                 "louvain12" = "Granus",
#                                 "louvain5" = "HSPCs",
#                                 "louvain1" = "HSPCs",
#                                 "louvain16" = "Basophils",
#                                 "louvain3" = "Basophils",
#                                 "louvain14" = "MiddleCells",
#                                 "louvain8" = "Bcells1",
#                                 "louvain10" = "Bcells2", 
#                                 "louvain9" = "Bcells3", 
#                                 "louvain17" = "NKs")

merge.louvs.active.hash <- list("louvain15" = "Eryths",
                                "louvain2" = "Eryths",
                                "louvain13" = "DCs",
                                "louvain6" = "DCs",
                                "louvain7" = "Granus",
                                "louvain4" = "Granus",
                                "louvain12" = "Granus",
                                "louvain5" = "HSPCs",
                                "louvain1" = "HSPCs",
                                "louvain16" = "Basophils",
                                "louvain3" = "Basophils",
                                "louvain14" = "MiddleCells",
                                "louvain8" = "Bcells",
                                "louvain10" = "Bcells", 
                                "louvain9" = "Bcells", 
                                "louvain17" = "NKs")


jhash <- hash::hash(merge.louvs.active.hash)

dat.umap.k4me1$cluster.annot.act <- sapply(dat.umap.k4me1$cluster, function(x) AssignHash(x = x, jhash = jhash, null.fill = NA))

ggplot(dat.umap.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.k4me1, aes(x = umap1, y = umap2, color = cluster.annot.act)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey95") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# show k9me3
dat.umap.k9me3$cluster.annot.act <- sapply(dat.umap.k9me3$louv.act, function(x) AssignHash(x = x, jhash = jhash, null.fill = NA))

ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = cluster.annot.act)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey95") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Filter out bad cells ----------------------------------------------------


varfilt.final <- 1.2
good.cells.dbl <- (dat.umap.merged.var %>% filter(stain == "dbl" & mark == "H3K9me3" & (cell.var.within.sum.norm >= varfilt.final | umap2 < -2)))$cell
good.cells <- c(good.cells.dbl, subset(dat.umap.merged.var, stain == "single")$cell)
                                            
dat.umap.merged.var.filt <- subset(dat.umap.merged.var, cell %in% good.cells)
                        
ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K9me3"), 
       aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")

ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K4me1"), 
       aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")


# Assign each repressed louvain to an active louvain ----------------------

dat.umap.merged.var.filt$cluster.act <- sapply(dat.umap.merged.var.filt$louv.act, function(x) AssignHash(x = x, jhash = jhash, null.fill = NA))

ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K9me3"), 
       aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")

ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K9me3"), 
       aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")

ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K4me1"), 
       aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")

ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K9me3"), 
       aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette, na.value = "grey95")



# Check plates ------------------------------------------------------------

dat.umap.merged.var.filt$plate <- sapply(dat.umap.merged.var.filt$cell, function(x) ClipLast(x, jsep = "_"))
dat.umap.merged.var.filt$experi <- sapply(dat.umap.merged.var.filt$plate, function(x) ClipLast(x, jsep = "-"))

ggplot(dat.umap.merged.var.filt %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~experi) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# assign k9me3 cluster to celltype by voting

jmark.tmp <- "H3K9me3"
inf.projs.tmp <- file.path(projmain, paste0("project_unmixed_", jmark.tmp, ".RData"))
load(inf.projs.tmp, v=T)

tm.result.tmp <- posterior(out.objs$out.lda)
tm.result.tmp <- AddTopicToTmResult(tm.result.tmp)

dat.umap.tmp <- DoUmapAndLouvain(tm.result.tmp$topics, jsettings = jsettings) %>%
  rowwise() %>%
  mutate(stain = "single")

tm.result.dbl.tmp <- AddTopicToTmResult(out.lda.predict)

dat.umap.dbl.tmp <- DoUmapAndLouvain(tm.result.dbl.tmp$topics, jsettings = jsettings) %>%
  rowwise() %>%
  mutate(stain = "dbl")


# Merge on one umap -------------------------------------------------------

umap.out.orig <- umap(tm.result.tmp$topics, config = jsettings)
umap.out.pred <- predict(umap.out.orig, data = tm.result.dbl.tmp$topics)

dat.umap.long.tmp1 <- data.frame(cell = rownames(umap.out.orig[["layout"]]), umap1 = umap.out.orig[["layout"]][, 1], umap2 = umap.out.orig[["layout"]][, 2], stain = "single", stringsAsFactors = FALSE)
dat.umap.long.tmp2 <- data.frame(cell = rownames(umap.out.pred), umap1 = umap.out.pred[, 1], umap2 = umap.out.pred[, 2], stain = "dbl", stringsAsFactors = FALSE)

dat.umap.long.tmp <- rbind(dat.umap.long.tmp1, dat.umap.long.tmp2)

topics.mat.tmp <- rbind(tm.result.tmp$topics, tm.result.dbl.tmp$topics)

dat.umap.long.tmp <- DoLouvain(topics.mat.tmp, jsettings, dat.umap.long.tmp)

dat.umap2.orig <- DoUmapAndLouvain(tm.result.tmp$topics, jsettings = jsettings)
dat.umap2.orig$louvorig <- dat.umap2.orig$louvain
cell2louvorig <- hash::hash(dat.umap2.orig$cell, dat.umap2.orig$louvorig)


ggplot(dat.umap2.orig, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("dat.umap2.orig original louvain numbers") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# remove bad cells from double stain
dat.umap.long.tmp <- subset(dat.umap.long.tmp, cell %in% good.cells)

dat.umap.long.tmp <- left_join(dat.umap.long.tmp, dat.umap2.orig)

dat.umap.long.tmp$louvorig <- sapply(dat.umap.long.tmp$cell, function(x) AssignHash(x, jhash = cell2louvorig, null.fill = NA))

jtmp <- dat.umap.long.tmp %>%
  filter(!is.na(louvorig)) %>%
  group_by(louvain, louvorig) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louvain) %>%
  mutate(ncells.frac = ncells / sum(ncells)) %>%
  filter(ncells == max(ncells))

louv2louvorig <- hash::hash(jtmp$louvain, jtmp$louvorig)

# table(subset(dat.umap.long.tmp, !is.na(louvorig), select = c(louvain, louvorig)))

ggplot(dat.umap.long.tmp, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~stain) + 
  ggtitle("Check louvain (ordering of clusters may differ because we have included projected cells)") + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# annotate cluster
dat.umap.long.tmp.merge <- left_join(dat.umap.long.tmp, subset(dat.umap.merged.var.filt, select = c(cell, cluster.act)))

ggplot(dat.umap.long.tmp.merge, aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~stain) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.long.tmp.merge, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~stain) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.umap.long.tmp.merge.sum.dbl <- dat.umap.long.tmp.merge %>%
  filter(stain == "dbl") %>%
  group_by(louvain, cluster.act) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louvain) %>%
  mutate(ncells.frac = ncells / sum(ncells))

# vote with weight
dat.umap.long.tmp.merge.sum.dbl.best <- dat.umap.long.tmp.merge.sum.dbl %>%
  group_by(louvain) %>%
  filter(ncells == max(ncells)) %>%
  ungroup() %>%
  filter(ncells > 5)

dat.umap.long.tmp.merge.sum.dbl.best$louvorig <- sapply(dat.umap.long.tmp.merge.sum.dbl.best$louvain, function(x) AssignHash(x = as.character(x), jhash = louv2louvorig, null.fill = NA))

louv.repress.hash <- hash::hash(paste("louvain", dat.umap.long.tmp.merge.sum.dbl.best$louvorig, sep = ""), dat.umap.long.tmp.merge.sum.dbl.best$cluster.act)

# annotate dat.umap.merged.var.filt single stains

dat.umap.merged.var.filt.annot <- dat.umap.merged.var.filt %>%
  rowwise() %>%
  mutate(cluster.act = ifelse(mark == "H3K4me1" & stain == "single", AssignHash(cluster, jhash, null.fill = NA), cluster.act),
         cluster.act = ifelse(mark == "H3K9me3" & stain == "single", AssignHash(cluster, louv.repress.hash, null.fill = NA), cluster.act))


ggplot(dat.umap.merged.var.filt.annot %>% filter(mark == "H3K4me1"), aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  ggtitle("Final K4me1 UMAP assigned by double staining")  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.merged.var.filt.annot %>% filter(mark == "H3K9me3"), aes(x = umap1, y = umap2, color = cluster.act)) + 
  geom_point() + 
  facet_wrap(~stain) + 
  ggtitle("Final K9me3 UMAP assigned by double staining") + 
  scale_color_manual(values = cbPalette, na.value = "grey95") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

dev.off()

# Save objects  -----------------------------------------------------------

save(dat.umap.merged.var.filt.annot, dat.umap.long.tmp.merge.sum.dbl.best, coords.dbl.annots, dat.umap.k4me1, dat.umap.k9me3, dat.umap.dbl.merge, good.cells, file = outfinal)
