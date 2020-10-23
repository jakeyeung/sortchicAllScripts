# Jake Yeung
# Date of Creation: 2020-10-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/compare_UMAPs_before_after_filtering.R
# Maybe I removed all HSPCs? 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf.k4me1.before <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl/cluster_tables_H3K4me1_BM_all_round2.txt"
inf.k9me3.before <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl/cluster_tables_H3K9me3_BM_all_round2.txt"

inf.k4me1.after <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl.cellfilt_binfilt/cluster_tables_H3K4me1_BM_all_round2.txt"
inf.k9me3.after <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl.cellfilt_binfilt/cluster_tables_H3K9me3_BM_all_round2.txt"

dat.k4me1.before <- fread(inf.k4me1.before)
dat.k4me1.after <- fread(inf.k4me1.after)

dat.k9me3.before <- fread(inf.k9me3.before)
dat.k9me3.after <- fread(inf.k9me3.after)

ggplot(dat.k9me3.before, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# cells.highlight <- subset(dat.k9me3.before, cluster == NA)$cell
cells.highlight <- subset(dat.k9me3.before, is.na(cluster))$cell

ggplot(dat.k9me3.after, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.k9me3.after %>% rowwise() %>% mutate(highlight = cell %in% cells.highlight), aes(x = umap1, y = umap2, color = highlight)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

