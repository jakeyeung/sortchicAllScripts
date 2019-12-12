# Jake Yeung
# Date of Creation: 2019-12-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/1-load_good_cells_filter_good_bins_other_winsize.R
# Automatically get good bins

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(scchicFuncs)
library(Matrix)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(tidyr)
library(GGally)
library(here)


# constants ---------------------------------------------------------------

outdir <- "~/data/scchic/quality_control_B6_other_winsizes"

blfile <- "~/data/scchic/databases/blacklist/mm10.blacklist.bed.gz"
assertthat::assert_that(file.exists(blfile))

dir.create(outdir)


# Load good cells  --------------------------------------------------------

jmark <- "H3K4me3"
TA.cutoff <- 0.5
count.cutoff <- 10^3
pcutoff.bin <- 0.95

inf.goodcells <- "/home/jyeung/data/scchic/quality_control_B6/JY_H3K4me3_bamnames.out"

cells.keep.dat <- read.table(inf.goodcells, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  dplyr::rename(cell = V1)

cells.keep <- cells.keep.dat$cell


# Load matrix -------------------------------------------------------------

indir.mat <- "/home/jyeung/data/from_cluster/scchic/ZellerRawDataB6_mergedAll.retag/countTables_otherWinSize"
infs.mat <- list.files(path = indir.mat, pattern = ".csv", full.names = TRUE, )

bsizes <- sapply(infs.mat, function(x) strsplit(basename(x), split = "\\.")[[1]][[4]], USE.NAMES = FALSE)
names(infs.mat) <- bsizes
names(bsizes) <- bsizes

mats <- lapply(infs.mat, function(x){
  ReadMatSlideWinFormat(x)
})

mats.filt <- lapply(mats, function(mat){
  mat[, cells.keep]
})


# Filter blacklist and correlated bins ------------------------------------

# rnames.all.common <- rownames(mats[[1]])

bl.dat <- LoadBlacklist(inf = blfile, asGR = FALSE)
bincor.dat <- read.table("~/data/scchic/quality_control_B6/B6_correlated_bins.bincutoff_0.95.2019-12-05.bed")
colnames(bincor.dat) <- c("seqnames", "start", "end")
bincor.dat$seqnames <- gsub("chrchr", "chr", bincor.dat$seqnames)
bl.dat.merge <- rbind(bl.dat, bincor.dat)
bl.gr <- makeGRangesFromDataFrame(bl.dat.merge)

FilterBlacklist <- function(rnames.all.common, bl.gr){
  rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all.common, GetChromo, add.chr=FALSE),
                                                   start = sapply(rnames.all.common, GetStart),
                                                   end = sapply(rnames.all.common, GetEnd)))
  overlaps <- findOverlaps(bl.gr, rnames.gr)
  indx <- seq(length(rnames.gr))
  bl.hits.i <- unique(subjectHits(overlaps))
  bl.hits.l <- !indx %in% bl.hits.i
  rnames.gr.filt <- rnames.gr[bl.hits.l]
  rnames.gr.badbins <- rnames.gr[!bl.hits.l]
  rnames.all.common.blfilt <- names(rnames.gr.filt)
  return(rnames.all.common.blfilt)
}

rnames.blfilt.lst <- lapply(mats.filt, function(mat){
  # get filtered blacklist
  rnames <- rownames(mat)
  return(FilterBlacklist(rnames, bl.gr))
})


# Do final bin filter on matrix -------------------------------------------

# filter out
mats.binfilt.lst <- lapply(bsizes, function(jname){
  return(mats.filt[[jname]][rnames.blfilt.lst[[jname]], ])
})

lapply(mats.binfilt.lst, dim)



# Measure sparsity --------------------------------------------------------

lapply(mats.binfilt.lst, function(mat){
  1 - Matrix::nnzero(mat) / length(mat)
})


# Write merged all --------------------------------------------------------

for (jwinsize in bsizes){
  mat <- mats.binfilt.lst[[jwinsize]]
  saveRDS(mat, file = file.path(outdir, paste0("B6BM_AllMerged_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".winsize_", jwinsize, ".rds")))
}

