# Jake Yeung
# Date of Creation: 2020-03-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/check_MARA_inputs_marksmerged.R
# Error in marks merged analysis check inputs

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


infN <- "/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged.txt"
N <- fread(infN)

infE <- "/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/count_mats_peaks_norm_merged/countmat_PZ_fromHiddenDomains_marks_merged.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11.txt"
E <- fread(infE)

rows.remove <- which(duplicated(N))

print(length(rows.remove))

N <- N[!duplicated(N)]

cnames.common <- intersect(N$Gene.ID, E$Gene.ID)

# check dimensions
Nmat <- subset(N, Gene.ID %in% cnames.common, select = -Gene.ID)
Emat <- subset(E, Gene.ID %in% cnames.common, select = -Gene.ID)
