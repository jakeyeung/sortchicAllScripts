t# Jake Yeung
# Date of Creation: 2022-01-05
# File: ~/projects/scchic/scripts/macbook_analysis_2021/public_data_analysis/check_Wu_2021_outputs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

inf.barcodes <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds"

dat.barcodes <- readRDS(inf.barcodes)
dat.barcodes.wu <- subset(dat.barcodes, grepl("Wu", jset))

ggplot(dat.barcodes.wu, aes(x = log10(nreads))) +
  geom_density()


inf.all <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/public_data/SRR12638101_1.sorted_dupcounts.bed"
inf.peaks <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/public_data/SRR12638101_1.sorted_dupcounts.hiddendomains_filt.bed"


# Set colnames ------------------------------------------------------------

cnames.all <- c("chromo", "start", "end", "allcounts", "bc", "umi")
cnames.peak <- c("chromo", "start", "end", "peakcounts", "bc", "umi")

# All ---------------------------------------------------------------------

dat.all <- fread(inf.all, col.names = cnames.all)

dat.all.sum <- dat.all %>%
  group_by(bc) %>%
  summarise(allcounts = sum(allcounts))


# Count peaks  ------------------------------------------------------------

dat.peaks <- fread(inf.peaks, col.names = cnames.peak)

dat.peaks.sum <- dat.peaks %>%
  group_by(bc) %>%
  summarise(peakcounts = sum(peakcounts))

dat.merged <- left_join(dat.all.sum, dat.peaks.sum)

minreads <- 0
ggplot(dat.merged %>% filter(allcounts >= minreads), aes(x = log10(allcounts))) +
  geom_density() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merged %>% filter(allcounts >= minreads), aes(x = log10(allcounts))) +
  geom_density() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merged %>% filter(allcounts >= minreads), aes(x = peakcounts / allcounts)) +
  geom_density() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merged %>% filter(allcounts >= minreads), aes(y = peakcounts / allcounts)) +
  geom_boxplot() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


