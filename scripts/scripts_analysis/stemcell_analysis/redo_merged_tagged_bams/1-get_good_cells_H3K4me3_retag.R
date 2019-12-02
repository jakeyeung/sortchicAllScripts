# Jake Yeung
# Date of Creation: 2019-11-30
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/redo_merged_tagged_bams/1-get_good_cells_H3K4me3_retag.R
# Retag and get a better estimate for cutoff

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)


GrepAndWriteMat <- function(mat.tmp, jgrp, jgrp.name, outf){
  cols.i <- grepl(jgrp, colnames(mat.tmp))
  mat.tmp.filt <- mat.tmp[, cols.i]
  assertthat::assert_that(ncol(mat.tmp.filt) > 0)
  print(jgrp.name)
  print(jgrp)
  print(dim(mat.tmp.filt))
  # write to output
  saveRDS(mat.tmp.filt, file = outf)
  return(mat.tmp.filt)
}


# Set cutoffs -------------------------------------------------------------

jmark <- "H3K4me3"
jmarks <- jmark; names(jmarks) <- jmark
cutoff.TA <- 0.5
cutoff.counts <- 1000
pcutoff <- paste0("TAcutoff_", cutoff.TA, ".", "cellsize_", cutoff.counts, ".", "oldbins")

outdir <- "/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged.abscutoff.dedupfixed"
dir.create(outdir)

# Read RZ -----------------------------------------------------------------

inf.rz <- "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells.retag/RZcounts/H3K4me3-BM_SC-merged.tagged.retagged.LH_counts.demuxbugfixed.csv"
inf.count <- "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells.retag/countTables/H3K4me3-BM_SC-merged.tagged.retagged.countTable.demuxbugfixed.csv"

empty.wells <- GetEmptyWells()

dat <- ReadLH.SummarizeTA(inf.rz, remove.nones = FALSE) %>%
  rowwise() %>% 
  mutate(TA.count = ifelse(is.na(TA.count), 0, TA.count), 
         TA.frac = ifelse(is.na(TA.frac), 0, TA.frac),
         experi = ClipLast(samp), 
         plate = strsplit(KeepLast(samp), "_")[[1]][[1]],
         experiplate = paste(experi, plate, sep = "_"), 
         cellindx = paste("cell", strsplit(KeepLast(samp), "_")[[1]][[2]], sep = ""),
         is.empty = cellindx %in% empty.wells,
         mark = GetMarkFromStr(experi))

ggplot(dat, aes(x = total.count, y = TA.frac)) + geom_point() + scale_x_log10() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat <- dat %>%
  rowwise() %>%
  mutate(good.cell = TA.frac > cutoff.TA & total.count > cutoff.counts)

ggplot(dat, aes(x = total.count, y = TA.frac, color = good.cell)) + geom_point() + scale_x_log10() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Use same bins as previous analysis  -------------------------------------

inf.old <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
load(inf.old, v=T)
mat.old <- count.dat$counts
rownames(mat.old) <- paste("chr", rownames(mat.old), sep = "")

bins.keep <- rownames(mat.old)

cells.keep <- subset(dat, good.cell)$samp

# load new mat ------------------------------------------------------------

mat <- ReadMatSlideWinFormat(inf.count)

rows.i <- which(rownames(mat) %in% bins.keep)
cols.i <- which(colnames(mat) %in% cells.keep)
mat.filt <- mat[rows.i, cols.i]


# Write all ---------------------------------------------------------------

mats.binfilt.lst <- list()
mats.binfilt.lst$H3K4me3 <- mat.filt

for (jmark in jmarks){
  print(jmark)
  mat <- mats.binfilt.lst[[jmark]]
  print(dim(mat))
  saveRDS(mat, file = file.path(outdir, paste0("B6BM_AllMerged_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds")))
}



# Write unenriched, lineage neg, stem cells, and merged -------------------

jmark <- "H3K4me3"

jgrp.bm <- paste0("^B6-13W1-BM-")
jgrp.linneg <- paste0("^PZ-Bl6-BM-Linneg-")
jgrp.sc <- paste0("^PZ-ChIC-B6BMSC-")
jgrp.vec <- c(jgrp.bm, jgrp.linneg, jgrp.sc)
jgrp.names <- c("Unenriched", "Linneg", "StemCells")
names(jgrp.vec) <- jgrp.names

for (jgrp.name in c("Unenriched", "Linneg", "StemCells")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}


# Write pairwise combinations  -----------------------------------------------------

# Unenriched X Linneg
jmarks.tmp <- c("H3K4me3")
# Unenriched, Linneg, Stemcell by itself
jgrp.bmXlinneg <- paste(jgrp.bm, jgrp.linneg, sep = "|")
jgrp.bmXlinneg.name <- "UnenrichedXLinneg"

for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXlinneg.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXlinneg, jgrp.bmXlinneg.name, outf)
}


# Unenriched X SC
jmarks.tmp <- c("H3K4me3")
jgrp.bmXsc <- paste(jgrp.bm, jgrp.sc, sep = "|")
jgrp.bmXsc.name <- "UnenrichedXStemCells"

# Unenriched, Linneg, Stemcell by itself
for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXsc.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXsc, jgrp.bmXsc.name, outf)
}

# Unenriched X SC
jmarks.tmp <- c("H3K4me3")

jgrp.linnegXsc <- paste(jgrp.linneg, jgrp.sc, sep = "|")
jgrp.linnegXsc.name <- "LinnegXStemCells"

# Unenriched, Linneg, Stemcell by itself
for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.linnegXsc.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.linnegXsc, jgrp.linnegXsc.name, outf)
}




