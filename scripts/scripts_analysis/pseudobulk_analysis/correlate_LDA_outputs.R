# Jake Yeung
# Date of Creation: 2019-04-01
# File: ~/projects/scchic/scripts/scripts_analysis/pseudobulk_analysis/correlate_LDA_outputs.R
# Do correlations on LDA outputs

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks.RData"

load(inf, v=T)


# Do correlations by topics?  ---------------------------------------------


# pick a reference cell: H3K4me1 erythryroblasts 
# topic 26?
jmark <- "H3K4me1"
jmark.compare <- "H3K9me3"

jtop <- "26"  # erythryoblast cells probably
jtop <- "13"  # 
topics.mat.h3k4me1 <- tm.result.lst[[jmark]]topics
eryth.weights <- sort(topics.mat.h3k4me1[, jtop], decreasing = TRUE)

# ref.cell <- "BM_H3K4me1_m2_rep1_cell141"
ref.cell <- names(eryth.weights)[[1]]

# check where this cell is on umap
ggplot(dat.merged.lst[[jmark]], aes(x = umap1, y = umap2, color = ifelse(cell == ref.cell, TRUE, FALSE))) + geom_point() 

ref.vec <- count.imputed.lst[[jmark]][, ref.cell]

dat.ref <- data.frame(logexprs.ref = ref.vec, bin = names(ref.vec))

# correlate with all cells in H3K9me3

compare.long <- gather(as.data.frame(count.imputed.lst[[jmark.compare]]) %>% mutate(bin = rownames(count.imputed.lst[[jmark.compare]])), key = cell, value = logexprs, -bin)
compare.long <- dplyr::left_join(compare.long, dat.ref, by = "bin")

compare.sum <- compare.long %>%
  group_by(cell) %>%
  summarise(correlation = cor(logexprs, logexprs.ref, method = "pearson", use = "complete.obs"))

ggplot(subset(compare.long, cell == "BM_H3K9me3_m1_rep1_cell1"), aes(x = logexprs, y = logexprs.ref)) + geom_point(alpha = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot UMAP, add color
dat.sub <- left_join(dat.merged.lst$H3K9me3 %>% filter(motif == "Ar"), compare.sum)

ggplot(dat.sub, aes(x = umap1, y = umap2, color = -correlation)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.sub %>% filter(correlation < -0.3), aes(x = umap1, y = umap2, color = correlation)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


