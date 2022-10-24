# Jake Yeung
# Date of Creation: 2022-05-21
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/37-replot_comparison_data.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf.rds <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/from_vivek/data_for_study_comparison.Rds"
dat.output <- readRDS(inf.rds)
top_bin_counts <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/from_vivek/top_regions_25p_count_fraction.Rds")

dat.output$top_bin_counts <- top_bin_counts

# Plot unique counts boxplot  ---------------------------------------------

dat.total <- dat.output$total_counts %>%
  bind_rows() %>%
  mutate(study = gsub("this_study", "zthis_study", study))


dat.gini <- dat.output$QC %>%
  bind_rows() %>%
  mutate(.id = gsub("this_study", "zthis_study", .id), 
         study = .id,
         type = sample)

dat.top <- dat.output$top_bin_counts %>%
  rowwise() %>%
  mutate(frac_in_top = X..i..)

# Unique reads ------------------------------------------------------------


# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.total, aes(x = study, y = uniq, fill = study)) + 
  geom_boxplot() + 
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~type, scales = "free_x", ncol = 1) + 
  theme_bw() + 
  xlab("") + 
  ylab("Unique fragments") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")


ggplot(dat.total, aes(x = study, y = uniq_total_ratio, fill = study)) + 
  geom_boxplot() + 
  # scale_y_log10(labels = function(x) format(x, scientific = FALSE)) + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~type, scales = "free_x", ncol = 1) + 
  theme_bw() + 
  xlab("") + 
  ylab("Fraction of Unique Fragments") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")


ggplot(dat.gini, aes(x = study, y = gini_coefficient, fill = study)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~type, scales = "free_x", ncol = 1) + 
  theme_bw() + 
  xlab("") + 
  ylab("Gini Coefficient") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")

# ggplot(dat.gini, aes(x = study, y = pct_counts_in_top_500_genes, fill = study)) + 
#   geom_boxplot() + 
#   scale_fill_manual(values = cbPalette) + 
#   facet_wrap(~type, scales = "free_x", ncol = 1) + 
#   theme_bw() + 
#   theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")

ggplot(dat.top, aes(x = study, y = frac_in_top, fill = study)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~type, scales = "free_x", ncol = 1) + 
  theme_bw() + 
  xlab("") + 
  ylab("Fraction of Fragments in Top 25% of Bins") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")






