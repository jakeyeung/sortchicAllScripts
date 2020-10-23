# Jake Yeung
# Date of Creation: 2020-08-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/check_TA_frac_with_without_spikeins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load before -------------------------------------------------------------



# Load after  -------------------------------------------------------------


indir.after <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/merged_bams/RZcounts.NewFilters.chr1only"
# indir.after <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/merged_bams/RZcounts.NewFilters.NoSpikeInChromo.manual_filt2"
indir.after2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/merged_bams/RZcounts.NewFilters"
# indir.after <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/merged_bams/RZcounts.NewFilters.SpikeInOnly"

infs.after <- list.files(path = indir.after, all.files = TRUE, pattern = "*.csv", full.names = TRUE)
names(infs.after) <- sapply(infs.after, function(x) strsplit( strsplit(basename(x), "\\.")[[1]][[1]], "-")[[1]][[3]])

infs.after2 <- list.files(path = indir.after2, all.files = TRUE, pattern = "*.csv", full.names = TRUE)
names(infs.after2) <- sapply(infs.after2, function(x) strsplit( strsplit(basename(x), "\\.")[[1]][[1]], "-")[[1]][[3]])

dat <- lapply(names(infs.after), function(jname){
  print(jname)
  inf <- infs.after[[jname]]
  dat.tmp <- ReadLH.SummarizeTA(inf)
  dat.tmp$mark <- jname
  dat.tmp$type <- "NoSpikeInChromo"
  return(dat.tmp)
}) %>%
  bind_rows()

dat2 <- lapply(names(infs.after2), function(jname){
  print(jname)
  inf <- infs.after2[[jname]]
  dat.tmp <- ReadLH.SummarizeTA(inf)
  dat.tmp$mark <- jname
  dat.tmp$type <- "All"
  return(dat.tmp)
}) %>%
  bind_rows()

print(range(dat$total.count))

ggplot(dat, aes(y = TA.frac, x = log10(total.count), color = mark)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  facet_wrap(~mark)

ggplot(dat2, aes(y = TA.frac, x = log10(total.count), color = mark)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  facet_wrap(~mark)

head(dat2)

cnames.keep <- intersect(dat$samp, dat2$samp)
plot(subset(dat, samp %in% cnames.keep)$total.count, subset(dat2, samp %in% cnames.keep)$total.count)
abline(a = 0, b = 1)
