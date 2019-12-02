# Jake Yeung
# Date of Creation: 2019-11-27
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/redo_merged_tagged_bams/1-get_good_cells.R
# Get good cells 

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

setwd(here())

# Constants ---------------------------------------------------------------

# same as paper 
# pcellsize <- 0.95
# pfrac <- 0.95
pcutoff.size <- 1000
pcutoff <- 0.5
pcutoff.bin <- 0.95

outdir <- "/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged.abscutoff.dedupbug"
dir.create(outdir)

# Read tables -------------------------------------------------------------

indir <- "/Users/yeung/data/scchic/from_cluster/ZellerRawDataB6_mergedAll/RZcounts"
indir.mat <- "/Users/yeung/data/scchic/from_cluster/ZellerRawDataB6_mergedAll/countTables"

jpatall <- paste0("*.csv")
jpatall.mats <- paste0("*.csv")
infs.rz.all <- list.files(indir, pattern = jpatall, full.names = TRUE)
infs.mats.all <- list.files(indir.mat, pattern = jpatall.mats, full.names = TRUE)

assertthat::assert_that(length(infs.rz.all) > 0)
assertthat::assert_that(length(infs.mats.all) > 0)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
infs.mats.all.name <- sapply(infs.mats.all, function(x) strsplit(basename(x), split = "-")[[1]][[1]])
names(infs.mats.all) <- infs.mats.all.name



# Plot counts -------------------------------------------------------------

empty.wells <- GetEmptyWells()
dats <- ReadLH.SummarizeTA(infs.rz.all, remove.nones = FALSE) %>%
  rowwise() %>% 
  mutate(TA.count = ifelse(is.na(TA.count), 0, TA.count), 
         TA.frac = ifelse(is.na(TA.frac), 0, TA.frac),
         experi = ClipLast(samp), 
         plate = strsplit(KeepLast(samp), "_")[[1]][[1]],
         experiplate = paste(experi, plate, sep = "_"), 
         cellindx = paste("cell", strsplit(KeepLast(samp), "_")[[1]][[2]], sep = ""),
         is.empty = cellindx %in% empty.wells,
         mark = GetMarkFromStr(experi))

dats <- dats %>%
  group_by(mark) %>%
  # ungroup() %>%
  mutate(good.cellsize = ifelse(total.count > pcutoff.size, TRUE, FALSE),
         good.frac = ifelse(TA.frac > 0.5, TRUE, FALSE),
         good.cell = as.logical(good.cellsize * good.frac * !is.empty))

m.qc <- ggplot(dats, aes(x = log10(total.count), y = TA.frac, shape = is.empty, size = is.empty, color = good.cell)) + geom_point() + 
  facet_wrap(experi ~ plate) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.qc)

# removed bad experiments
min.good.cells <- 150
good.cells.by.plate <- dats %>%
  group_by(experiplate) %>%
  summarise(n.goodcells = length(which(good.cell))) %>%
  arrange(n.goodcells)

bad.plates <- subset(good.cells.by.plate, n.goodcells < min.good.cells)

print(paste("removing bad plates:", bad.plates$experiplate))

dats <- subset(dats, !experiplate %in% bad.plates)

cells.keep <- subset(dats, good.cell)$samp  # filtered out bad plates already 


# Load mats ---------------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
mats <- lapply(infs.mats.all, function(inf){
  jmat <- ReadMatSlideWinFormat(inf)
  # filter bad chromosomes
  
  rows.i <- sapply(rownames(jmat), function(x) strsplit(x, ":")[[1]][[1]] %in% jchromos)
  cols.i <- which(colnames(jmat) %in% cells.keep)
  if (length(cols.i) > 0){
    return(jmat[rows.i, cols.i])
  } else{
    return(NULL)
  }
})

# remove NULLs in list
mats <- purrr::compact(mats)


# merge by mark
mats.merge <- lapply(jmarks, function(jmark){
  mats.sublst.i <- which(names(mats) == jmark)
  mats.filt <- mats[mats.sublst.i]
  
  # get temporary common rows (AFTER removing bad cells)
  rows.common.temp <- Reduce(intersect, lapply(mats.filt, rownames))
  mats.filt.filt <- lapply(mats.filt, function(x){
    rows.i <- which(rownames(x) %in% rows.common.temp)
    return(x[rows.i, ])
  }) %>%
    do.call(cbind, .)
})

print(lapply(mats.merge, dim))

rnames.all.common <- Reduce(intersect, lapply(mats, function(mat) rownames(mat)))


# Get bin means -----------------------------------------------------------

bin.lst <- lapply(jmarks, function(jmark){
  Matrix::rowMeans(mats.merge[[jmark]])
})


# Filter blacklist --------------------------------------------------------

rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all.common, GetChromo, add.chr=FALSE),
                                                 start = sapply(rnames.all.common, GetStart),
                                                 end = sapply(rnames.all.common, GetEnd)))
bl.gr <- LoadBlacklist()

overlaps <- findOverlaps(bl.gr, rnames.gr)

indx <- seq(length(rnames.gr))
bl.hits.i <- unique(subjectHits(overlaps))
bl.hits.l <- !indx %in% bl.hits.i
rnames.gr.filt <- rnames.gr[bl.hits.l]
rnames.gr.badbins <- rnames.gr[!bl.hits.l]


# Use filtered blacklist --------------------------------------------------

rnames.all.common.blfilt <- names(rnames.gr.filt)



# Find correlations across marks ------------------------------------------

bin.long.nofilt <- lapply(jmarks, function(jmark){
  x <- bin.lst[[jmark]]
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
}) %>%
  bind_rows()
bin.wide.nofilt <- spread(bin.long.nofilt, key = mark, value = bin.mean)

ggpairs(log10(bin.wide.nofilt %>% dplyr::filter(coord %in% names(rnames.gr.badbins)) %>% dplyr::select(-coord)), 
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Blacklist bin signal")

# filter out blacklist before doing correlations
bin.long <- lapply(jmarks, function(jmark){
  x <- bin.lst[[jmark]]
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
  dat <- subset(dat, coord %in% rnames.all.common.blfilt)
}) %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(bin.rank = rank(bin.mean) / length(bin.mean))

bin.sum <- bin.long %>%
  group_by(coord) %>%
  summarise(bin.rank.mean = mean(bin.rank),
            bin.rank.min = min(bin.rank),
            bin.rank.max = max(bin.rank))


bin.wide <- spread(bin.long %>% dplyr::select(-bin.rank), key = mark, value = bin.mean)

ggpairs(log10(bin.wide %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Before removing correlated bins")

# remove bad bins?
bins.high <- subset(bin.sum, bin.rank.mean >= pcutoff.bin)$coord
bins.low <- subset(bin.sum, bin.rank.mean <= (1-pcutoff.bin))$coord

bins.correlated <- c(bins.high, bins.low)
print(length(bins.correlated))

ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("After removing correlated bins")

# check any systematic biases on the offset
bsize <- 100000
bin.sum$bin <- as.numeric(sapply(bin.sum$coord, GetStart))
bin.sum$offset <- sapply(bin.sum$bin, function(x) x %% bsize)

ggplot(bin.sum, aes(x = bin.rank.mean, group = offset)) + geom_histogram() + facet_wrap(~offset, ncol = 1) +
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



plot(table(sapply(bins.high, GetChromo)))
plot(table(sapply(bins.low, GetChromo)))

rnames.blfilt.corrfilt <- (bin.wide %>% dplyr::filter(!coord %in% bins.correlated))$coord


unique(sapply(rnames.blfilt.corrfilt, function(x) strsplit(x, ":")[[1]][[1]]))


# Write correlated bins to output -----------------------------------------

# bins.correlated.dat <- data.frame(bin = bins.correlated)
bins.correlated.dat <- data.frame(chromo = paste0("chr", sapply(bins.correlated, GetChromo)),
                                  start = sapply(bins.correlated, GetStart),
                                  end = sapply(bins.correlated, GetEnd))

outbed <- file.path(outdir, paste0("B6_correlated_bins.", Sys.Date(), ".bed"))
fwrite(bins.correlated.dat, file = outbed, col.names = FALSE, sep = "\t")

for (jmark in jmarks){
  cells.keep.bymark <- data.frame(cell = paste0(grep(pattern = jmark, cells.keep, value = TRUE), ".sorted.bam"))
  fwrite(cells.keep.bymark, file = file.path(outdir, paste0("JY_", jmark, "_bamnames.out")), col.names = FALSE, sep = "\t")
}


# write pdf
pdf(file = file.path(outdir, paste0("qc_plots_B6_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".pdf")))
print(m.qc)

ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("After removing correlated bins")

ggpairs(log10(bin.wide %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Before removing correlated bins")
dev.off()


# Do final bin filter on matrix -------------------------------------------

# filter out
mats.binfilt.lst <- lapply(mats.merge, function(mat){
  return(mat[rnames.blfilt.corrfilt, ])
})

lapply(mats.binfilt.lst, dim)


# Measure sparsity --------------------------------------------------------

lapply(mats.binfilt.lst, function(mat){
  1 - Matrix::nnzero(mat) / length(mat)
})


# Write merged all --------------------------------------------------------



for (jmark in jmarks){
  print(jmark)
  mat <- mats.binfilt.lst[[jmark]]
  print(dim(mat))
  saveRDS(mat, file = file.path(outdir, paste0("B6BM_AllMerged_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds")))
}



# Write single unenriched -------------------------------------------------


# H3K4me1 different from rest?
# H3K4me1 something weird, skip for now 
# unique(sapply(colnames(mats.binfilt.lst$H3K4me1), ClipLast))

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

# generl grep strings no mark
# jgrp.bm <- paste0("^B6-13W1-BM-", jmark)
# jgrp.linneg <- paste0("^PZ-Bl6-BM-Linneg-", jmark)
# jgrp.sc <- paste0("^PZ-ChIC-B6BMSC-", jmark)
# jgrp.vec <- c(jgrp.bm, jgrp.linneg, jgrp.sc)
# jgrp.names <- c("Unenriched", "Linneg", "StemCells")
# names(jgrp.vec) <- jgrp.names

jgrp.bm <- paste0("^B6-13W1-BM-")
jgrp.linneg <- paste0("^PZ-Bl6-BM-Linneg-")
jgrp.sc <- paste0("^PZ-ChIC-B6BMSC-")
jgrp.vec <- c(jgrp.bm, jgrp.linneg, jgrp.sc)
jgrp.names <- c("Unenriched", "Linneg", "StemCells")
names(jgrp.vec) <- jgrp.names


jmark <- "H3K4me3"
for (jgrp.name in c("Unenriched", "Linneg", "StemCells")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}

jmark <- "H3K27me3"
for (jgrp.name in c("Unenriched", "Linneg", "StemCells")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}

jmark <- "H3K9me3"
for (jgrp.name in c("Unenriched", "Linneg")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}


# Write pairwise combinations  -----------------------------------------------------

# Unenriched X Linneg
jmarks.tmp <- c("H3K4me3", "H3K27me3", "H3K9me3")
# Unenriched, Linneg, Stemcell by itself
jgrp.bmXlinneg <- paste(jgrp.bm, jgrp.linneg, sep = "|")
jgrp.bmXlinneg.name <- "UnenrichedXLinneg"

for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXlinneg.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXlinneg, jgrp.bmXlinneg.name, outf)
}


# Unenriched X SC
jmarks.tmp <- c("H3K4me3", "H3K27me3")
jgrp.bmXsc <- paste(jgrp.bm, jgrp.sc, sep = "|")
jgrp.bmXsc.name <- "UnenrichedXStemCells"

# Unenriched, Linneg, Stemcell by itself
for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXsc.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXsc, jgrp.bmXsc.name, outf)
}

# Unenriched X SC
jmarks.tmp <- c("H3K4me3", "H3K27me3")

jgrp.linnegXsc <- paste(jgrp.linneg, jgrp.sc, sep = "|")
jgrp.linnegXsc.name <- "LinnegXStemCells"

# Unenriched, Linneg, Stemcell by itself
for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.linnegXsc.name, "_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.linnegXsc, jgrp.linnegXsc.name, outf)
}



