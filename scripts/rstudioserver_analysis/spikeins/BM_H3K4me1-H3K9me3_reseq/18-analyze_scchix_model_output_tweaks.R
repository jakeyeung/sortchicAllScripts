# Jake Yeung
# Date of Creation: 2021-01-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/18-analyze_scchix_model_output.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load output -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

wmin <- 0.1
wmax <- 0.9
pcount <- 0

# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.4.wmax_0.6.RData")
# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.wmax_0.7.RData")
# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.wmax_0.7.pseudocount_0.RData")

infrdata <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_", wmin, ".wmax_", wmax, ".pseudocount_", pcount, ".RData"))

assertthat::assert_that(file.exists(infrdata))


# "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.4.wmax_0.6.RData"

load(infrdata, v=T)

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

library(ggforce)
m.grid <- ggplot(coords.dbl.annots, aes(x = louv.act, y = louv.repress, louv.act, color = lnprob)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")

print(m.grid)



# Check UMAPs  ------------------------------------------------------------


# load UMAP 
inf.umap <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/lda_and_clusters/ClusterAnnot.lda_and_datmerged.k4_k9_dynamic_regions.H3K4me1xH3K9me3.2021-01-30.RData"
load(inf.umap, v=T)

dat.umap.dbl <- dat.merge

dat.umap.dbl.merge <- left_join(dat.umap.dbl, coords.dbl.annots) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

m.check.repress <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

m.check.act <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = louv.act)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

m.check.w <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = w)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_viridis_c()

JFuncs::multiplot(m.check.repress, m.check.act, cols = 2)

print(m.check.w)

ggplot(dat.umap.dbl.merge %>% filter(w < 0.8 & w > 0.2), aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

ggplot(dat.umap.dbl.merge %>% filter(w < 0.7), aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 


# Check total cuts --------------------------------------------------------

cuts.total <- data.frame(cell = colnames(count.mat), cuts = colSums(count.mat), stringsAsFactors = FALSE)

dat.umap.dbl.merge2 <- left_join(dat.umap.dbl.merge, cuts.total)

ggplot(dat.umap.dbl.merge2, aes(x = umap1, y = umap2, color = log2(cuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.dbl.merge2, aes(x = umap1, y = umap2, color = w)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(density(log2(dat.umap.dbl.merge2$cuts)))


jcutoff <- 2^8
# lncutoff <- -0.00001
# lncutoff <- -0.01
lncutoff <- -10^-4
wcutoff <- 0.7

m.check.repress2  <- ggplot(dat.umap.dbl.merge2 %>% filter(cuts > jcutoff & lnprob > lncutoff), aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

m.check.act2 <- ggplot(dat.umap.dbl.merge2 %>% filter(cuts > jcutoff & lnprob > lncutoff), aes(x = umap1, y  = umap2, color = louv.act)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

JFuncs::multiplot(m.check.repress2, m.check.act2, cols = 2)

ggplot(dat.umap.dbl.merge2 %>% filter(cuts > jcutoff & lnprob > lncutoff & w < wcutoff), aes(x = umap1, y  = umap2, color = louv.repress)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85")  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.umap.dbl.merge2 %>% mutate(louv.repress = ifelse(cuts > jcutoff & lnprob > lncutoff & w < wcutoff, louv.repress, NA)), aes(x = umap1, y  = umap2, color = louv.repress)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85")  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")




m.grid2 <- ggplot(dat.umap.dbl.merge2 %>% filter(cuts > jcutoff & lnprob > lncutoff), aes(x = louv.act, y = louv.repress, louv.act, color = lnprob)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")

print(m.grid2)

# 
# 
# # Check fits on erythroblasts  --------------------------------------------
# 
# jsub <- subset(dat.umap.dbl.merge, louv.act == "Eryths" & louv.repress != "Eryths")
# cells.eryth <- jsub$cell
# 
# jcell <- cells.eryth[[4]]
# fits.out[[jcell]]
# # subset(coords.dbl, )
# 
# 
