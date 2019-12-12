# Jake Yeung
# Date of Creation: 2019-12-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/1-filter_mat_good_cells_good_bins.R
# Good cells, good bins

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

# Functions ---------------------------------------------------------------

GrepAndWriteMat <- function(mat.tmp, jgrp, jgrp.name, outf){
  cols.i <- grepl(jgrp, colnames(mat.tmp))
  mat.tmp.filt <- mat.tmp[, cols.i]
  if (ncol(mat.tmp.filt) == 0){
    warning("Empty matrix after filtering, returning NA")
    return(NA)
  }
  print(jgrp.name)
  print(jgrp)
  print(dim(mat.tmp.filt))
  print(1 - Matrix::nnzero(mat.tmp.filt) / length(mat.tmp.filt))
  # write to output
  saveRDS(mat.tmp.filt, file = outf)
  return(mat.tmp.filt)
}


# Cutoffs -----------------------------------------------------------------

jchromos <- paste("chr", c(seq(25)), sep = "")  # Zebrafish 

# outdir <- "/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged.abscutoff.dedupbug.allmarks"
outdir <- "/home/jyeung/data/scchic/quality_control_ZF"
dir.create(outdir)
blfile <- "~/data/scchic/databases/blacklist/mm10.blacklist.bed.gz"
assertthat::assert_that(file.exists(blfile))

dir.create(outdir)
TA.cutoff <- 0.5
count.cutoff <- 500
pcutoff.bin <- 0.95

# Load counts -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/data/from_cluster/scchic/ZellerRawDataZF_all.retag"
dir.rz <- file.path(indir, "RZcounts")
dir.counts <- file.path(indir, "countTables")


infs.rz <- list.files(path = dir.rz, pattern = "*.csv", full.names = TRUE)
infs.counts <- list.files(path = dir.counts, pattern = "*.csv", full.names = TRUE)

assertthat::assert_that(length(infs.rz) > 0)
assertthat::assert_that(length(infs.counts) > 0)


names(infs.rz) <- sapply(infs.rz, function(x) strsplit(basename(x), split = "-")[[1]][[1]])
names(infs.counts) <- sapply(infs.counts, function(x) strsplit(basename(x), split = "-")[[1]][[1]])

dat.rz <- lapply(infs.rz, function(inf){
  jtmp <- ReadLH.SummarizeTA(inf, remove.nones = FALSE) %>%
    rowwise() %>%
    mutate(mark = GetMarkFromStr(samp),
           experi = ClipLast(samp))
  return(jtmp)
}) %>%
  bind_rows()

# remove -G1 from experi name
dat.rz$experi <- gsub("-G1$", "", dat.rz$experi)

dat.rz$mark <- factor(as.character(dat.rz$mark), levels = jmarks)


# Get good cells ----------------------------------------------------------

empty.wells <- GetEmptyWells()

dat.rz <- dat.rz %>%
  rowwise() %>%
  mutate(good.cell = TA.frac >= TA.cutoff & total.count >= count.cutoff,
         cellindx = paste("cell", sapply(samp, function(x){
           jsplit <- strsplit(x, "_")[[1]]
           return(jsplit[[length(jsplit)]])
         }), sep = ""), 
         is.empty = cellindx %in% empty.wells)


cells.keep <- subset(dat.rz, good.cell)$samp

print(paste("Number of cells before:", nrow(dat.rz)))
print(paste("Number of cells after:", length(cells.keep)))

# Plot them ---------------------------------------------------------------

# plot them
m.all <- ggplot(dat.rz, aes(x = total.count, y = TA.frac, color = good.cell, shape = is.empty, size = is.empty)) + geom_point() + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10()

# plot each mark separately, then show experiments
lapply(jmarks, function(jmark){
  m <- ggplot(dat.rz %>% filter(mark == jmark), aes(x = total.count, y = TA.frac, color = good.cell, shape = is.empty, size = is.empty)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi) + 
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
  print(m)
})

print(unique(dat.rz$experi))

# Load mats ---------------------------------------------------------------

mats <- lapply(infs.counts, function(inf){
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

rnames.all.common <- Reduce(intersect, lapply(mats, function(mat) rownames(mat)))



# Get bin means -----------------------------------------------------------

bin.lst <- lapply(jmarks, function(jmark){
  Matrix::rowMeans(mats[[jmark]])
})


# Filter blacklist --------------------------------------------------------

rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all.common, GetChromo, add.chr=FALSE),
                                                 start = sapply(rnames.all.common, GetStart),
                                                 end = sapply(rnames.all.common, GetEnd)))
bl.gr <- makeGRangesFromDataFrame(data.frame(seqnames = "chrXXX", start = 123, end = 567))
# # bl.gr <- LoadBlacklist(inf = blfile)
# bl.gr <- makeGRangesFromDataFrame(list())  # no blacklist for zebrafish yet
# 
overlaps <- findOverlaps(rnames.gr, bl.gr)
# 
indx <- seq(length(rnames.gr))
bl.hits.i <- unique(subjectHits(overlaps))
bl.hits.l <- !indx %in% bl.hits.i
rnames.gr.filt <- rnames.gr[bl.hits.l]
rnames.gr.badbins <- rnames.gr[!bl.hits.l]

# rnames.gr.filt <- rnames.gr


# Use filtered blacklist --------------------------------------------------

rnames.all.common.blfilt <- names(rnames.gr.filt)



# Plot bins ---------------------------------------------------------------

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

m.pairs.clean <- ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
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

print(length(rnames.blfilt.corrfilt))



# Write correlated bins to output -----------------------------------------

# plot each mark separately, then show experiments
lapply(jmarks, function(jmark){
  m <- ggplot(dat.rz %>% filter(mark == jmark), aes(x = total.count, y = TA.frac)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi) + 
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
  print(m)
})

ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
        lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("After removing correlated bins")

ggpairs(log10(bin.wide %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
  ggtitle("Before removing correlated bins")
dev.off()


# Do final bin filter on matrix -------------------------------------------

# filter out
mats.binfilt.lst <- lapply(mats, function(mat){
  return(mat[rnames.blfilt.corrfilt, ])  # cells already filtered
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
  print(1 - Matrix::nnzero(mat) / length(mat))
  saveRDS(mat, file = file.path(outdir, paste0("ZF_AllMerged_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds")))
}


# Write singles -----------------------------------------------------------


jgrp.bm <- paste0("^PZ-ChIC-ZFWKM-")
jgrp.sc <- paste0("^PZ-ChIC-ZFWKMCD41plus-")
jgrp.vec <- c(jgrp.bm, jgrp.sc)
jgrp.names <- c("Unenriched", "StemCells")
names(jgrp.vec) <- jgrp.names

# jmark <- "H3K4me3"

for (jmark in jmarks){
  print(jmark)
  for (jgrp.name in c("Unenriched", "StemCells")){
    print(jgrp.name)
    jgrp <- jgrp.vec[[jgrp.name]]
    outf <- file.path(outdir, paste0("ZF_", jgrp.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
    GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
  }
}

# Write pairwise combinations  -----------------------------------------------------

jmarks.tmp <- c("H3K4me3", "H3K4me3", "H3K9me3")  # no H3K27me3 yet
jgrp.bmXsc <- paste(jgrp.bm, jgrp.sc, sep = "|")
jgrp.bmXsc.name <- "UnenrichedXStemCells"

for (jmark in jmarks.tmp){
  print(jmark)
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("ZF_", jgrp.bmXsc.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXsc, jgrp.bmXsc.name, outf)
}


