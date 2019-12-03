# Jake Yeung
# Date of Creation: 2019-12-03
# File: ~/projects/scchic/scripts/scripts_analysis/all_merged_analysis/1-get_good_cells_merged_dedupbugfixed.R
# Debugging


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


# Functions ---------------------------------------------------------------

GrepAndWriteMat <- function(mat.tmp, jgrp, jgrp.name, outf){
  cols.i <- grepl(jgrp, colnames(mat.tmp))
  mat.tmp.filt <- mat.tmp[, cols.i]
  assertthat::assert_that(ncol(mat.tmp.filt) > 0)
  print(jgrp.name)
  print(jgrp)
  print(dim(mat.tmp.filt))
  print(1 - Matrix::nnzero(mat.tmp.filt) / length(mat.tmp.filt))
  # write to output
  saveRDS(mat.tmp.filt, file = outf)
  return(mat.tmp.filt)
}


# Cutoffs -----------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/quality_control_PZ_Bl6-BM-Linneg-Stemcells_AllMerged.abscutoff.dedupbug.allmarks"
dir.create(outdir)
TA.cutoff <- 0.5
count.cutoff <- 10^3
pcutoff.bin <- 0.95

# Load counts -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

dir.rz <- "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells.retag.demuxbugfixed/RZcounts"
dir.counts <- "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells.retag.demuxbugfixed/countTables"

infs.rz <- list.files(path = dir.rz, pattern = "*.csv", full.names = TRUE)
infs.counts <- list.files(path = dir.counts, pattern = "*.csv", full.names = TRUE)

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

cells.keep <- dat.rz$samp

# Plot them ---------------------------------------------------------------

# plot them
m.all <- ggplot(dat.rz, aes(x = total.count, y = TA.frac, color = good.cell, shape = is.empty, size = is.empty)) + geom_point() + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10()

# plot each mark separately, then show experiments
lapply(jmarks, function(jmark){
  m <- ggplot(dat.rz %>% filter(mark == jmark), aes(x = total.count, y = TA.frac)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi) + 
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
  print(m)
})


# Load mats ---------------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
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
bl.gr <- LoadBlacklist()

overlaps <- findOverlaps(bl.gr, rnames.gr)

indx <- seq(length(rnames.gr))
bl.hits.i <- unique(subjectHits(overlaps))
bl.hits.l <- !indx %in% bl.hits.i
rnames.gr.filt <- rnames.gr[bl.hits.l]
rnames.gr.badbins <- rnames.gr[!bl.hits.l]


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

# bins.correlated.dat <- data.frame(bin = bins.correlated)
bins.correlated.dat <- data.frame(chromo = paste0("chr", sapply(bins.correlated, GetChromo)),
                                  start = sapply(bins.correlated, GetStart),
                                  end = sapply(bins.correlated, GetEnd))

outbed <- file.path(outdir, paste0("B6_correlated_bins.bincutoff_", pcutoff, ".", Sys.Date(), ".bed"))
fwrite(bins.correlated.dat, file = outbed, col.names = FALSE, sep = "\t")

for (jmark in jmarks){
  cells.keep.bymark <- data.frame(cell = paste0(grep(pattern = jmark, cells.keep, value = TRUE), ".sorted.bam"))
  fwrite(cells.keep.bymark, file = file.path(outdir, paste0("JY_", jmark, "_bamnames.out")), col.names = FALSE, sep = "\t")
}


# write pdf
pdf(file = file.path(outdir, paste0("qc_plots_B6.", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".pdf")))
print(m.all)

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
  saveRDS(mat, file = file.path(outdir, paste0("B6BM_AllMerged_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds")))
}


# Write singles -----------------------------------------------------------

# H3K4me1 is different because we used a new antibody? 
jgrp.bm.h3k4me1 <- paste0("^PZ-ChIC-Bl6-BM-H3K4me1-Index")
jgrp.sc.h3k4me1 <- paste0("^PZ-ChIC-Bl6-BM-stem-cells-H3K4me1-Index")

jgrp.bm <- paste0("^B6-13W1-BM-")
jgrp.linneg <- paste0("^PZ-Bl6-BM-Linneg-")
jgrp.sc <- paste0("^PZ-ChIC-B6BMSC-")
jgrp.vec <- c(jgrp.bm, jgrp.linneg, jgrp.sc)
jgrp.names <- c("Unenriched", "Linneg", "StemCells")
names(jgrp.vec) <- jgrp.names

jgrp.vec.h3k4me1 <- c(jgrp.bm.h3k4me1, jgrp.sc.h3k4me1)
names(jgrp.vec.h3k4me1) <- c("Unenriched", "StemCells")

jmark <- "H3K4me1"
for (jgrp.name in c("Unenriched", "StemCells")){
  jgrp <- jgrp.vec.h3k4me1[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}

jmark <- "H3K4me3"
for (jgrp.name in c("Unenriched", "Linneg", "StemCells")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}

jmark <- "H3K27me3"
for (jgrp.name in c("Unenriched", "Linneg", "StemCells")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}

jmark <- "H3K9me3"
for (jgrp.name in c("Unenriched", "Linneg")){
  jgrp <- jgrp.vec[[jgrp.name]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp, jgrp.name, outf)
}


# Write pairwise combinations  -----------------------------------------------------

# Unenriched X Linneg
jmarks.tmp <- c("H3K4me3", "H3K27me3", "H3K9me3")  # no H3K4me1 yet
# Unenriched, Linneg, Stemcell by itself
jgrp.bmXlinneg <- paste(jgrp.bm, jgrp.linneg, sep = "|")
jgrp.bmXlinneg.name <- "UnenrichedXLinneg"

for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXlinneg.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXlinneg, jgrp.bmXlinneg.name, outf)
}


# Unenriched X SC: for H3K4me3 and H3K27me3
jmarks.tmp <- c( "H3K4me3", "H3K27me3")
jgrp.bmXsc <- paste(jgrp.bm, jgrp.sc, sep = "|")
jgrp.bmXsc.name <- "UnenrichedXStemCells"

for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXsc.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXsc, jgrp.bmXsc.name, outf)
}

# Unenriched X SC: for H3K4me1
jmarks.tmp <- c( "H3K4me1")
jgrp.bmXsc <- paste(jgrp.bm.h3k4me1, jgrp.sc.h3k4me1, sep = "|")
jgrp.bmXsc.name <- "UnenrichedXStemCells"

for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.bmXsc.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.bmXsc, jgrp.bmXsc.name, outf)
}


# Linneg X SC (H3K4me1 excluded)
jmarks.tmp <- c("H3K4me3", "H3K27me3")
jgrp.linnegXsc <- paste(jgrp.linneg, jgrp.sc, sep = "|")
jgrp.linnegXsc.name <- "LinnegXStemCells"

# Unenriched, Linneg, Stemcell by itself
for (jmark in jmarks.tmp){
  mat.tmp <- mats.binfilt.lst[[jmark]]
  outf <- file.path(outdir, paste0("B6BM_", jgrp.linnegXsc.name, "_", jmark, ".TAcutoff_", TA.cutoff, ".countscutoff_", count.cutoff, ".binfilt_cellfilt.", Sys.Date(), ".rds"))
  GrepAndWriteMat(mats.binfilt.lst[[jmark]], jgrp.linnegXsc, jgrp.linnegXsc.name, outf)
}




