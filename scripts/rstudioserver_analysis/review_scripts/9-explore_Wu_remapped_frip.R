# Jake Yeung
# Date of Creation: 2021-06-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/9-explore_Wu_remapped_frip.R
# Remap Wu, call peaks. Check outputs

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Check before  -----------------------------------------------------------


inf.barcodes <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds"
dat.barcodes <- readRDS(inf.barcodes)
dat.barcodes.wu <- subset(dat.barcodes, grepl("Wu", jset))

ggplot(dat.barcodes.wu, aes(x = log10(nreads))) + 
  geom_density() 



# Count all  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

cnames.all <- c("chromo", "start", "end", "allcounts", "bc", "umi")
cnames.peak <- c("chromo", "start", "end", "peakcounts", "bc", "umi")

inf.all <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/counts_output_r1_notrim_bowtie2/SRR12638101_1.sorted_dupcounts.bed")
dat.all <- fread(inf.all, col.names = cnames.all)

dat.all.sum <- dat.all %>%
  group_by(bc) %>%
  summarise(allcounts = sum(allcounts))

# Count peaks  ------------------------------------------------------------

inf.peaks <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/peak_filtered_beds_again/SRR12638101_1.sorted_dupcounts.hiddendomains_filt.bed")
dat.peaks <- fread(inf.peaks, col.names = cnames.peak)

dat.peaks.sum <- dat.peaks %>%
  group_by(bc) %>%
  summarise(peakcounts = sum(peakcounts))

dat.merged <- left_join(dat.all.sum, dat.peaks.sum)

minreads <- 30
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



