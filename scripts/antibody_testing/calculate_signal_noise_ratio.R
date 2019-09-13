# Jake Yeung
# Date of Creation: 2019-08-01
# File: ~/projects/scchic/scripts/antibody_testing/calculate_signal_noise_ratio.R
# Calculate signal to noise ratio across antibody concentrations and batches




library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)


# Load data  --------------------------------------------------------------

# indir <- "/Users/yeung/data/scchic/bigwigs/antibody_tests"

# inf.new <- "/Users/yeung/data/scchic/bigwigs/antibody_tests/PZ_K4me1_four_results.txt"
inf.new <- "/Users/yeung/data/scchic/bigwigs/antibody_tests/PZ_K4me1_four_results_100kb.txt"
# inf.orig <- "/Users/yeung/data/scchic/bigwigs/antibody_tests/H3K4me1_ChIP_100kb_original_4.txt"
# inf.orig <- "/Users/yeung/data/scchic/bigwigs/antibody_tests/PZ-K4me1-75-1-4_100kb.txt"
inf.orig <- "/Users/yeung/data/scchic/bigwigs/antibody_tests/PZ-K4me1-75-1-4_100kb_100kbbins.txt"

dat.new <- fread(inf.new); colnames(dat.new)[[1]] <- "chromo"; colnames(dat.new) <- gsub("'", "", colnames(dat.new)); colnames(dat.new) <- gsub("_", "-", colnames(dat.new))

dat.orig <- fread(inf.orig); colnames(dat.orig)[[1]] <- "chromo"; colnames(dat.orig) <- gsub("'", "", colnames(dat.orig)); colnames(dat.orig) <- gsub("_", "-", colnames(dat.orig))

chromos.new <- unique(dat.new$chromo)
chromos.old <- unique(dat.orig$chromo)

chromos.common <- intersect(chromos.new, chromos.old)

dat.new <- subset(dat.new, chromo %in% chromos.common)
dat.orig <- subset(dat.orig, chromo %in% chromos.common)

dat.merge <- left_join(dat.new, dat.orig, by = c("chromo", "start", "end"))
# remove rows containing NaNs
good.rows <- apply(dat.merge, MARGIN = 1, FUN = function(jrow) all(jrow != "NaN" & !is.na(jrow)))
dat.merge <- dat.merge[good.rows, ]

dat.merge.long <- dat.merge %>%
  tidyr::gather(key = "sample", value = "count", c(-chromo, -start, -end))

dat.merge.long$batch <- sapply(dat.merge.long$sample, function(s) strsplit(s, "-")[[1]][[3]])
dat.merge.long$conc <- sapply(dat.merge.long$sample, function(s) strsplit(s, "-")[[1]][[5]])

# plot(density(dat.merge.long$count))
plot(density(log10(dat.merge.long$count)))

m.hist <- ggplot(dat.merge.long, aes(x = count)) + geom_histogram(bins = 500) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(conc ~ batch) + ggtitle("Distribution of counts across 100 kb bins")  + 
  scale_x_log10()
print(m.hist)

