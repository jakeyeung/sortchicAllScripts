# Jake Yeung
# Date of Creation: 2019-12-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/0-qcllines_in_bam_fastq.R
# Check we didnt hav some bug

# HVG-scChIC-Bl6-intestine-k9me3-index1-20-09-19 is bad run, it was sorting and tagging before it was done mapping!!


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf.sum <- "/home/jyeung/hpc/intestinal_scchic/raw_data/count_summaries.txt"

dat <- fread(inf.sum, header = FALSE) %>%
  rowwise() %>%
  mutate(fname = strsplit(V1, "\\.")[[1]][[1]],
         ftype = strsplit(V1, "\\.")[[1]][[2]],
         counts.corrected = ifelse(endsWith(ftype, "bam"), V2 / 2, V2 / 4))

dat.sum <- dat %>%
  group_by(fname) %>%
  filter(endsWith(ftype, "bam")) %>%
  summarise(counts.diff = diff(range(counts.corrected))) %>%
  arrange(desc(counts.diff))
  
dat.sum.fastq <- dat %>%
  group_by(fname) %>%
  filter(!endsWith(ftype, "bam")) %>%
  summarise(counts.diff = diff(range(counts.corrected)) / counts.corrected[[2]]) %>%
  arrange(desc(counts.diff))

dat.sum.bamfastq <- dat %>%
  group_by(fname) %>%
  filter(ftype %in% c("bam", "trim")) %>%
  summarise(counts.diff = diff(range(counts.corrected)) / counts.corrected[[2]]) %>%
  arrange(desc(counts.diff))
  
# BAD RUN
