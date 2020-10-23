# Jake Yeung
# Date of Creation: 2020-09-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/celltyping_with_chix/1-get_bins_from_cleaned_data.again.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load data  --------------------------------------------------------------

# inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
inf.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_H3K27me3_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
load(inf.lda, v=T)
mat.AllMerged <- count.mat
rnames.AllMerged <- rownames(mat.AllMerged)
# bins.keep <- rownames(count.mat.AllMerged)


# Load bins from Maria  ---------------------------------------------------

inf.dbl <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/2020-06-04_redo_count_tables/countTables.unfixed.50000_50000/all_BM_K27m3_200119.mq_40.bsize_50000.step_50000.csv.gz")
assertthat::assert_that(file.exists(inf.dbl))
mat <- ReadMatSlideWinFormat(inf.dbl, as.sparse = TRUE)


# Get good cells for ChIX  ---------------------------------------------------------

inf.lda.dbl <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_BM_dblmarks_deeper_2020-01-19.remove_bad_plates/lda_outputs.qc_removed_bad_plates.mark_K27m3.countcutoff_500.TAcutoff_0.5.redo.K-30.binarize.FALSE/ldaOut.qc_removed_bad_plates.mark_K27m3.countcutoff_500.TAcutoff_0.5.redo.K-30.Robj")
assertthat::assert_that(file.exists(inf.lda.dbl))

load(inf.lda.dbl, v=T)

count.mat.dbl.filt <- count.mat

good.cells <- colnames(count.mat.dbl.filt)


# Get common bins ---------------------------------------------------------

mat.filt <- mat[, good.cells]
rnames.ChIX <- rownames(mat.filt)


# Get common bins and rewrite ---------------------------------------------

bins.common <- intersect(rnames.AllMerged, rnames.ChIX)



# Merge everything and rerun LDA  -----------------------------------------

mat.merged <- cbind(mat.AllMerged[bins.common, ], mat.filt[bins.common, ])

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix"


# write to output

saveRDS(mat.merged, file = file.path(outdir, "H3K27me3_AllMerged_and_ChIX_merged.rds"))



# Also try projecting -----------------------------------------------------


# need to pad zeros 

rows.to.add <- setdiff(rnames.AllMerged, bins.common)

mat.to.add <- Matrix(data = 0, nrow = length(rows.to.add), ncol = ncol(mat.filt), dimnames = list(rows.to.add, colnames(mat.filt)))

mat.filt.padded <- rbind(mat.to.add, mat.filt[bins.common, ])
mat.filt.padded <- mat.filt.padded[rnames.AllMerged, ]  # reorder

saveRDS(mat.filt.padded, file = file.path(outdir, "H3K27me3_ChIX_padded.rds"))

