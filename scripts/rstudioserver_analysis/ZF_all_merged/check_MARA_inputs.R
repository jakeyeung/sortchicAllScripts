# Jake Yeung
# Date of Creation: 2020-08-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/check_MARA_inputs.R
# 

rm(list=ls())

library(scchicFuncs)

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/H3K4me1/mara_input/count_mats_peaks_norm/ldaOut.H3K4me1.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.keepNbins_250.txt"
inf2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/H3K4me1/mara_input/count_mats_peaks_norm/ldaOut.H3K4me1.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.txt"
inf.site <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/H3K4me1/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.H3K4me1.txt"

exprs <- read.table.handlerows(inf)
exprs2 <- read.table.handlerows(inf2)
site <- read.table.handlerows(inf.site)

site2 <- read.table(inf.site, header = TRUE)

exprs.sub <- read.table(inf, header = TRUE, nrows = 200)
site.sub <- read.table(inf.site, header = TRUE, nrows=200)

rnames.exprs <- make.names(exprs.sub[, 1], unique = TRUE)
rnames.site <- make.names(site.sub[, 1], unique = TRUE)
