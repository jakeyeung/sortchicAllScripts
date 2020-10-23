# Jake Yeung
# Date of Creation: 2020-08-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/remove_large_peaks_write_clean_N_and_E.BM.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(mixtools)

# Explore peaks  ----------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  print(jmark)
  
  inf.mat <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".txt"))
  outf.mat <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".Cleaned.txt"))
  outf.mat.colscaled <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".Cleaned.ColumnScaled.txt"))
  outf.mat.rowscaled <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".Cleaned.RowScaled.txt"))
  # outpdf <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".Cleaned.pdf"))
  
  # assertthat::assert_that(!file.exists(outf.mat))
  # assertthat::assert_that(!file.exists(outf.mat.colscaled))
  # assertthat::assert_that(!file.exists(outf.mat.rowscaled))
  # assertthat::assert_that(!file.exists(outpdf))
  
  
  mat <- read.table.handlerows(inf.mat)
  
  inf.exprs <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/count_mats_peaks_from_GLMPCA/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.txt"))
  # assertthat::assert_that(!file.exists(outf.exprs))
  
  
  E <- read.table.handlerows(inf.exprs)
  
  # Get common rows  --------------------------------------------------------
  
  rows.common <- intersect(rownames(mat), rownames(E))
  
  rows.keep <- rows.common
  
  
  
  # Write new sitecount mat and expr mat ------------------------------------
  
  Nmat <- mat[rows.keep, ]
  N <- data.frame(Gene.ID = rows.keep, Nmat, stringsAsFactors = FALSE)
  N.colscaled <- data.frame(Gene.ID = rownames(Nmat), as.matrix(scale(Nmat, center = FALSE, scale = TRUE)), stringsAsFactors = FALSE)
  N.rowscaled <- data.frame(Gene.ID = rownames(Nmat), as.matrix(t(scale(t(Nmat), center = FALSE, scale = TRUE))), stringsAsFactors = FALSE)
  
  Enew <- data.frame(Gene.ID = rows.keep, E[rows.keep, ], stringsAsFactors = FALSE)
  
  # Write new  --------------------------------------------------------------
  
  fwrite(N, file = outf.mat, sep = "\t")
  fwrite(N.colscaled, file = outf.mat.colscaled, sep = "\t")
  fwrite(N.rowscaled, file = outf.mat.rowscaled, sep = "\t")
  # fwrite(Enew, file = outf.exprs, sep = "\t")
  
}




