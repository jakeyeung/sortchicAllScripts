# Jake Yeung
# Date of Creation: 2019-06-10
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/analysis_checks/calculate_median_counts_Grosselin_et_al_2019.R
# Claculatet median counts per cell from Nature genetics paper

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)


# Load data tables  -------------------------------------------------------

indir <- "/Users/yeung/data/scchic/public_data/scChIPseq_data_Grosselin_et_al_2019"

# (fpaths <- list.files(indir, pattern = "*scChIP*", full.names = TRUE))
fpaths <- list.files(indir, pattern = "*scChIP.*K4me3.*", full.names = TRUE)
assertthat::assert_that(length(fpaths) > 0)

dats <- lapply(fpaths, function(x){
  dat <- fread(x)
  cname1 <- colnames(dat)[[1]]
  rownames(dat) <- dat[[cname1]]
  dat[[cname1]] <- NULL
  return(as.matrix(dat))
})

dats.size.k4me3 <- lapply(dats, colSums)

plot(density(log10(unlist(dats.size.k4me3))))

# median
print(median(unlist(dats.size.k4me3)))
mean(unlist(dats.size.k4me3))


fpaths <- list.files(indir, pattern = "*scChIP.*K27me3_hg38.txt", full.names = TRUE)
assertthat::assert_that(length(fpaths) > 0)

dats <- lapply(fpaths, function(x){
  dat <- fread(x)
  cname1 <- colnames(dat)[[1]]
  rownames(dat) <- dat[[cname1]]
  dat[[cname1]] <- NULL
  return(as.matrix(dat))
})

dats.size.k27me3 <- lapply(dats, colSums)

plot(density(log10(unlist(dats.size.k27me3))))

# median
print(median(unlist(dats.size.k27me3)))
print(mean(unlist(dats.size.k27me3)))
