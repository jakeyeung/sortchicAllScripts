# Jake Yeung
# Date of Creation: 2021-01-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/2-downstream_DE_analysis_HSPC_peaks_all_marks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load outputs ------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again.peaks"


out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir, paste0("poisson_fit_HSPCpeaks.", jmark, ".2021-01-15.newannot2.witherrors.RData"))
  load(inf.tmp, v=T)
  # fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  # fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-29.newannot2.witherrors.MoreBins.", jsuffix, ".RData")
  # fname <- paste0("poisson_fit_bins.", jmark, ".2020-12-29.newannot2.witherrors.MoreBins.HSPCs_vs_nonHSPCs.spikeins.RData")
  # inf.fits <- file.path(indir.fits, fname)
  load(inf.tmp, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, pvals.long = pvals.long))
})

params.lst <- lapply(out.lst, function(out){
  out$params.long
})

pvals.lst <- lapply(out.lst, function(out){
  out$pvals.long
})

ggplot(params.lst$H3K27me3, aes(x = estimate, y = -log10(pval.param), color = param)) + 
  facet_wrap(~param) + 
  geom_point() 

ggplot(params.lst$H3K9me3, aes(x = estimate, y = -log10(pval.param), color = param)) + 
  facet_wrap(~param) + 
  geom_point() 

# Check Gbe1  -------------------------------------------------------------





# Check Tead1 ?  ----------------------------------------------------------





