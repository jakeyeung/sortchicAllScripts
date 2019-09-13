# Jake Yeung
# Date of Creation: 2019-06-18
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/3-LDA_downstream.R
# LDA downstream

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)

data(AssociatedPress)

train <- AssociatedPress[1:100, ]
test <- AssociatedPress[101:150, ]

train.mat <- Matrix::Matrix(as.matrix(train), sparse=TRUE)
test.mat <- Matrix::Matrix(as.matrix(test), sparse=TRUE)

train.lda <- LDA(train.mat, 5)
(train.topics <- topics(train.lda))

test.topics <- posterior(train.lda, test.mat)


# Load original H3K4me3 LDA and project new data --------------------------

inf <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.RData"
load(inf, v=T)

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)

train.lda <- out.objs$out.lda


# Load new sparse matrix --------------------------------------------------

inf.test <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-Linneg/PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.2019-06-17.RData"
assertthat::assert_that(file.exists(inf.test))
load(inf.test, v=T)

test.topics <- posterior(train.lda, t(as.matrix(count.dat$counts)))


# Plot new samples --------------------------------------------------------

te



