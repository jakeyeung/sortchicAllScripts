# Jake Yeung
# Date of Creation: 2018-12-21
# File: ~/projects/scChiC/scripts/scripts_analysis/quality_controls/explore_peak_sizes.R
# After merging macs outputs by a distance using bedtools, check the size of peaks are proper 

library(dplyr)
library(data.table)
# inf <- "/private/tmp/macs2out/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_50000.broadPeak"
# inf <- "/private/tmp/macs2out/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_10000.broadPeak"
inf <- "/private/tmp/macs2out/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_25000.broadPeak"
# inf <- "/private/tmp/macs2out/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_100000.broadPeak"

cnames <- c("chromo", "start", "end", "name", "blank", "score", "strand", "rgb", "pval2", "qval")
dat <- data.table::fread(inf, col.names = cnames) %>% 
  mutate(dists = end - start)

ggplot(dat, aes(x = dists)) + geom_histogram() + scale_x_log10()
ggplot(dat, aes(x = dists)) + geom_histogram()

print(dim(dat))