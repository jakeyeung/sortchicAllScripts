# Jake Yeung
# Date of Creation: 2020-12-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/4-check_mean_effect.R
# Check mean effect (HSPC effect)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load  -------------------------------------------------------------------

indir.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3"

params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster")
})

pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

# pvalcutoff <- 10^-10



# Check skewness across marks  --------------------------------------------

params.sum.long <- lapply(jmarks, function(jmark){
  params.sum <- params.lst[[jmark]] %>%
    group_by(bin) %>%
    summarise(estimate.sum = sum(estimate) / length(estimate),
              se.sum = sqrt(sum(se ^ 2) / length(se))) %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

m <- ggplot(params.sum.long %>% filter(abs(estimate.sum) < 5), aes(x = estimate.sum, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)
