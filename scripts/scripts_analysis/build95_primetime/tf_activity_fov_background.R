# Jake Yeung
# Date of Creation: 2019-04-11
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/tf_activity_fov_background.R
# Calculate fov background and also the real fovs

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)


# Functions ---------------------------------------------------------------

GetPvalFOV <- function(fovs.permute, fovs.real = NULL, jprob = 0.9, show.plot = FALSE, return.pval.only = FALSE){
  jmark <- fovs.permute$mark[[1]]
  if (!is.null(fovs.real)){
    fov.real <- subset(fovs.real, mark == jmark)$fov
  } else {
    assertthat::assert_that(length(unique(fovs.permute$fov.real)) == 1)
    fov.real <- fovs.permute$fov.real[[1]]
  }
  
  fovs.sub <- fovs.permute %>%
    # filter(mark == jmark) %>%
    group_by(mark) %>%
    arrange(fov) %>%
    ungroup() %>%
    mutate(fov.cumsum = cumsum(fov),
           indx = seq(length(fov)),
           frac.less.than = indx / length(indx),
           frac.more.than = 1 - frac.less.than,
           log10.frac.more.than = log10(frac.more.than))
  
  
  fovs.subsub <- fovs.sub %>% filter(fov > quantile(fov, probs = jprob) & frac.more.than > 0)
  jfit <- lm(formula = log10.frac.more.than ~ fov, data = fovs.subsub)
  
  log10pval <- predict(jfit, newdata = data.frame(fov = fov.real))
  
  xpred <- seq(min(fovs.subsub$fov), max(fov.real, fovs.subsub$fov), length.out = 100)
  ypred <- predict(jfit, newdata = data.frame(fov = xpred))
  pred.dat <- data.frame(log10.frac.more.than = ypred, fov = xpred)
  
  if (show.plot){
    m <- ggplot(fovs.subsub, aes(x = fov, y = log10.frac.more.than)) + 
      geom_point() + theme_bw() +
      geom_vline(xintercept = fov.real, linetype = "dashed") + 
      expand_limits(y = ceiling(log10pval)) + 
      geom_line(mapping = aes(x = fov, y = log10.frac.more.than), data = pred.dat) + ggtitle(jmark)
    print(m)
  }
  if (!return.pval.only){
    return(list(real.dat = fovs.subsub, pred.dat = pred.dat, fit = jfit, log10pval = log10pval))
  } else {
    return(data.frame(log10pval = log10pval))
  }
}


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



