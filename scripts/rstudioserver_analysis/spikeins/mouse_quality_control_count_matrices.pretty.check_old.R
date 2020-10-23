# Jake Yeung
# Date of Creation: 2020-08-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse_quality_control_count_matrices.pretty.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(irlba)


inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40/countTablesAndRZr1only.NewFilters/H3K4me3-BM_AllMerged.merged_by_clusters_with_NAs.count_table.binsize_50000.csv"
mat <- ReadMatSlideWinFormat(inf)
plot(density(log10(colSums(mat))))

dat.old <- data.frame(cuts.total = colSums(mat), cell = colnames(mat), experi = "old", stringsAsFactors = FALSE)

inf.new <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.binsize_50000.csv"
mat.new <- ReadMatSlideWinFormat(inf.new)
plot(density(log10(colSums(mat.new))))

dat.new <- data.frame(cuts.total = colSums(mat.new), cell = colnames(mat.new), experi = "WithSpikeins", stringsAsFactors = FALSE)

dat.merge <- bind_rows(dat.old, dat.new)

ggplot(dat.merge, aes(x = cuts.total, fill = experi)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() 

ggplot(dat.merge, aes(x = log2(cuts.total), fill = experi)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = seq(15)) + 
  geom_vline(xintercept = 12)


