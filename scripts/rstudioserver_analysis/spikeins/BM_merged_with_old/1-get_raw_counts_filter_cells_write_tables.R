# Jake Yeung
# Date of Creation: 2020-10-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/1-get_raw_counts_filter_cells_write_tables.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load raw counts old -----------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
indir.old <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs")
indir.new <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.blfix/varfilt")

mats.old <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("count_mat_cbfilt_maxcountsfilt.all.2020-06-16.", jmark, ".rds")
  inf <- file.path(indir.old, fname)
  mat <- readRDS(inf)
})




# Load raw counts new from rds --------------------------------------------

mats.new <- lapply(jmarks, function(jmark){
  print(jmark)
  # fname <- paste0("count_mat_cbfilt_maxcountsfilt.all.2020-06-16.", jmark, ".rds")
  fname <- paste0("count_mat.", jmark, ".filt_0.15_0.95_counts_and_l2r.blfix.varfilt.rds")
  inf <- file.path(indir.new, fname)
  mat <- readRDS(inf)
})



# Match rows and merge ----------------------------------------------------

mats.merge <- lapply(jmarks, function(jmark){
  m1 <- mats.old[[jmark]]
  m2 <- mats.new[[jmark]]
  rows.common <- intersect(rownames(m1), rownames(m2))
  m.merge <- cbind(m1[rows.common, ], m2[rows.common, ])
})



# Write mats --------------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM"

lapply(jmarks, function(jmark){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_old_merged_with_new.", jmark, ".rds"))
  print(dim(mats.merge[[jmark]]))
  saveRDS(mats.merge[[jmark]], outf)
})

