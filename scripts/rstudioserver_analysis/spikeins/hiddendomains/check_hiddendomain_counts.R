# Jake Yeung
# Date of Creation: 2020-12-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/hiddendomains/check_hiddendomain_counts.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


hubprefix <- "/home/jyeung/hub_oudenaarden"


clusters <- c("Basophils", "Bcells", "DCs", "Eryths", "Granulocytes", "HSPCs", "NKs", "pDCs")
names(clusters) <- clusters

mlength <- "50000"
indir <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/hiddendomains_outputs/hiddendomains_outputs_minlength_", mlength))
assertthat::assert_that(dir.exists(indir))
infs <- lapply(clusters, function(clst){
  dname <- paste0("PZ-BM-H3K27me3-", clst, "-merged.", mlength, ".cutoff")
  inf <- file.path(indir, dname, paste0("PZ-BM-H3K27me3-", clst, "-merged.", mlength, ".cutoff_treatment_bins.txt"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

dat.lst <- lapply(clusters, function(clst){
  dat.tmp <- fread(infs[[clst]])
  dat.tmp$cluster <- clst
  return(dat.tmp)
})

dat.merged <- dat.lst %>% 
  bind_rows()

ggplot(dat.merged, aes(x = log10(count), fill = cluster)) + 
  facet_wrap(~cluster) + 
  geom_density() + 
  geom_vline(xintercept = c(log10(60), log10(10))) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

# plot(density(log2(dat.lst$Eryths$count)))
# plot(density(log2(dat.lst$pDCs$count)))
# plot(density(log2(dat.lst$NKs$count)))
# plot(density(log2(dat.lst$Granulocytes$count)))
# 

