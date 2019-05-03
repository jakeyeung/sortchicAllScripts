# Jake Yeung
# Date of Creation: 2019-04-11
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/tf_activity_fov_background.R
# Calculate fov background and also the real fovs

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)


# Functions ---------------------------------------------------------------

source("scripts/Rfunctions/BackgroundPermutationScripts.R")


# Load fovs ---------------------------------------------------------------

inf.fovs <- "/Users/yeung/data/scchic/from_cluster/fov_permute_summary/FOVs_summary.parsed.txt.gz"
fovs.permute <- read.table(gzfile(inf.fovs), header = FALSE, sep = "\t")
colnames(fovs.permute) <- c("fov", "seed", "mark")

marks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
indir.fovs.real <- c("/Users/yeung/data/scchic/from_cluster/mara_analysis_build95.cells_from_bin_analysis_FINAL/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50",
                   "/Users/yeung/data/scchic/from_cluster/mara_analysis_build95.cells_from_bin_analysis_FINAL/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me3.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50",
                   "/Users/yeung/data/scchic/from_cluster/mara_analysis_build95.cells_from_bin_analysis_FINAL/hiddenDomains_cellmin_0-cellmax_9999999-binarize_FALSE-BM_H3K27me3.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50",
                   "/Users/yeung/data/scchic/from_cluster/mara_analysis_build95.cells_from_bin_analysis_FINAL/hiddenDomains_cellmin_0-cellmax_9999999-binarize_FALSE-BM_H3K9me3.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50")
names(indir.fovs.real) <- marks

fovs.lst <- lapply(indir.fovs.real, function(indir){
  indir2 <- list.dirs(indir, full.names = TRUE)[[2]]
  inf <- list.files(indir2, pattern = "FOV", full.names = TRUE)
  assertthat::assert_that(file.exists(inf))
  dat <- unlist(data.table::fread(inf), use.names = FALSE)
  return(dat)
})
fovs <- data.frame(fov.real = unlist(fovs.lst), mark = marks)


# Check background model  -------------------------------------------------

jmark <- "H3K4me1"
jprob <- 0.9

jtitle <- jmark

fovs.permute <- left_join(fovs.permute, fovs)
pvals.long <- fovs.permute %>%
  group_by(mark) %>%
  do(GetPvalFOV(., fovs.real = NULL, jprob = 0.9, show.plot = TRUE, return.pval.only = TRUE))
# ggplot(fovs.sub, aes(x = fov, y = log10.frac.more.than)) + geom_point()  + facet_wrap(~mark)



