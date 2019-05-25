# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/3b-get_count_sums_from_LDA_input.R
# Calculate count sums from LDA non-binarized input 

library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")

ReadLDAGetCountMat <- function(jmark, inf){
  print(paste("Reading from", inf))
  out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = "auto")
  return(out.objs$count.mat)
}


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt"

infs <- lapply(jmarks, function(jmark, d){
  return(list.files(d, pattern = paste0("_", jmark, "_.*K-30_40_50.Robj"), full.names = TRUE))
}, indir)

count.mat.lst <- lapply(jmarks, function(jmark) ReadLDAGetCountMat(jmark, infs[[jmark]]))


# Get cell sums  ----------------------------------------------------------

count.sum <- lapply(jmarks, function(jmark){
  x <- count.mat.lst[[jmark]]
  csums <- Matrix::colSums(x)
  csums.long <- data.frame(cell = names(csums), cellsum = csums, stringsAsFactors = FALSE, mark = jmark)
}) %>%
  bind_rows()


# Write to file  ----------------------------------------------------------

robjsdir <- "/Users/yeung/data/scchic/robjs/B6_objs"

save(count.sum, file = file.path(robjsdir, "cell_sums_from_LDA_input.RData"))


