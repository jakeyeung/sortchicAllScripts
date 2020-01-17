# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/check_MT_reads.R
# Check MT reads

library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf <- "/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_otherWinSize_NoSliding/H3K4me1-BM_SC-merged.tagged.bsize_50000.step_50000.countTable.demuxbugfixed.csv"
mat <- ReadMatSlideWinFormat(inf, as.sparse = TRUE)

mat.sub <- mat[grepl("chrM", rownames(mat)), ]
chromos <- unique(sapply(rownames(mat), function(x) strsplit(x, ":")[[1]][[1]]))
