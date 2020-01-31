# Jake Yeung
# Date of Creation: 2020-01-30
# File: ~/projects/scchic/scripts/macbook_analysis/make_bed_from_count_mat/create_bed_from_count_mat.R
# description

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(scchicFuncs)
library(JFuncs)

# Load count mat ----------------------------------------------------------

inf <- "/Users/yeung/data/scchic/from_cluster/LDA_outputs/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"
load(inf, v=T)

# Write bedfiles to output ------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# create bedfile
dat <- data.frame(chromo = sapply(rownames(count.mat), GetChromo),
                  start = sapply(rownames(count.mat), GetStart),
                  end = sapply(rownames(count.mat), GetEnd), 
                  stringsAsFactors = FALSE)
# filter relevant chromosomes 
dat <- subset(dat, chromo %in% jchromos)

dat <- dat %>%
  mutate(chromo = gsub("chr", "", chromo)) %>%
  arrange(chromo, as.numeric(start), as.numeric(end))

outf <- "/Users/yeung/data/scchic/beds_from_count_mat/BM_H3K4me1_all_bin_regions_50000_chromofilt.txt"

fwrite(dat, file = outf, sep = "\t", col.names = FALSE)