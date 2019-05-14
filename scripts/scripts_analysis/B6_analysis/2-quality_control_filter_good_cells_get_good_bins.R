# Jake Yeung
# Date of Creation: 2019-05-07
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/1-quality_control_on_bins.R
# Quality control of bins afterwards


library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(tidyr)
library(GGally)

LoadMats <- function(infs, cells.keep){
  sparse.mats <- lapply(infs, function(inf) readRDS(inf))
  # filter cells
  # sparse.mats <- lapply(sparse.mats, function(mat){
  #   cnames <- colnames(mat)
  #   cnames.keep.i <- which(cnames %in% cells.keep)
  #   return(mat[, cnames.keep.i])
  # })
  rnames <- lapply(sparse.mats, rownames)
  rnames.common <- Reduce(intersect, rnames)
  sparse.mats <- lapply(sparse.mats, function(x){
    return(x[rnames.common, ])
  })
  mat.merge <- do.call(cbind, sparse.mats)
  cells.keep.i <- which(colnames(mat.merge) %in% cells.keep)
  mat.merge <- mat.merge[, cells.keep.i]
  return(mat.merge)
}

LoadBlacklist <- function(inf = "data/blacklists/mm10.blacklist.bed.gz", asGR = TRUE){
  dat <- fread(inf, col.names = c("seqnames", "start", "end"))
  if (asGR){
    dat <- makeGRangesFromDataFrame(dat)
  }
  return(dat)
}

source("scripts/Rfunctions/Aux.R")

# Constants ---------------------------------------------------------------


# pcutoff <- 0.99
pcutoff <- 0.95
pcutoff.cell <- 0.05

prefix <- ".windowbug2"
outdir <- paste0("~/data/scchic/count_mat_binfilt_cellfilt_for_LDA", prefix)
dir.create(outdir)

# Load  -------------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

indir <- paste0("/Users/yeung/data/scchic/from_cluster/count_mat_from_bam", prefix)

infs.h3k4me1 <- list.files(indir, pattern = ".*H3K4me1.*.rds", full.names = TRUE)
infs.h3k4me3 <- list.files(indir, pattern = ".*H3K4me3.*.rds", full.names = TRUE)
infs.h3k27me3 <- list.files(indir, pattern = ".*H3K27me3.*.rds", full.names = TRUE)
infs.h3k9me3 <- list.files(indir, pattern = ".*H3K9me3.*.rds", full.names = TRUE)
infs.all <- list.files(indir, pattern = "*.rds", full.names = TRUE)
infs.bymark <- list(infs.h3k4me1, infs.h3k4me3, infs.h3k27me3, infs.h3k9me3)
names(infs.bymark) <- jmarks

inf.cells <- paste0("/Users/yeung/data/scchic/quality_control", prefix, "/good_cells.pcounts.", pcutoff.cell, ".pfrac.", pcutoff.cell, ".txt")
assertthat::assert_that(file.exists(inf.cells))
cells <- fread(inf.cells)
cells.keep <- cells$cell

bin.means.h3k4me1 <- Matrix::rowMeans(LoadMats(infs.h3k4me1, cells.keep))
bin.means.h3k4me3 <- Matrix::rowMeans(LoadMats(infs.h3k4me3, cells.keep))
bin.means.h3k27me3 <- Matrix::rowMeans(LoadMats(infs.h3k27me3, cells.keep))
bin.means.h3k9me3 <- Matrix::rowMeans(LoadMats(infs.h3k9me3, cells.keep))

# mats.test <- LoadMats(infs.h3k4me1, cells.keep)

# make matrix then ggpairs

bin.lst <- list(bin.means.h3k4me1, bin.means.h3k4me3, bin.means.h3k27me3, bin.means.h3k9me3)
names(bin.lst) <- jmarks

rnames.all <- lapply(bin.lst, function(x) names(x))

rnames.all.common <- Reduce(intersect, rnames.all)


# Filter blacklist --------------------------------------------------------

rnames.gr <- makeGRangesFromDataFrame(data.frame(seqnames = sapply(rnames.all.common, GetChromo, add.chr=TRUE),
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
  # dat <- subset(dat, coord %in% rnames.all.common.blfilt)
}) %>%
  bind_rows()
bin.wide.nofilt <- spread(bin.long.nofilt, key = mark, value = bin.mean)
ggpairs(log10(bin.wide.nofilt %>% filter(coord %in% names(rnames.gr.badbins)) %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() + 
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
bins.high <- subset(bin.sum, bin.rank.mean >= pcutoff)$coord
bins.low <- subset(bin.sum, bin.rank.mean <= (1-pcutoff))$coord

# bins.high <- subset(bin.sum, bin.rank.min >= pcutoff)$coord
# bins.low <- subset(bin.sum, bin.rank.max <= (1 - pcutoff))$coord

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


# Write correlated bins to output -----------------------------------------

# bins.correlated.dat <- data.frame(bin = bins.correlated)
bins.correlated.dat <- data.frame(chromo = paste0("chr", sapply(bins.correlated, GetChromo)),
                                 start = sapply(bins.correlated, GetStart),
                                 end = sapply(bins.correlated, GetEnd))

outbed <- file.path(outdir, paste0("B6_correlated_bins.", Sys.Date(), ".bed"))
fwrite(bins.correlated.dat, file = outbed, col.names = FALSE, sep = "\t")

for (jmark in jmarks){
  cells.keep.bymark <- data.frame(cell = paste0(grep(pattern = jmark, cells.keep, value = TRUE), ".sorted.bam"))
  fwrite(cells.keep.bymark, file = paste0("/Users/yeung/data/scchic/tables/bamlist_for_peak_analysis_build95_B6/JY_", jmark, "_bamnames.out"), col.names = FALSE, sep = "\t")
}

# write cells


# Do final bin filter on matrix -------------------------------------------

mats.lst <- lapply(infs.bymark, LoadMats, cells.keep)

# filter out
mats.binfilt.lst <- lapply(mats.lst, function(mat){
  return(mat[rnames.blfilt.corrfilt, ])
})

lapply(mats.binfilt.lst, dim)


# Measure sparsity --------------------------------------------------------

lapply(mats.binfilt.lst, function(mat){
  1 - Matrix::nnzero(mat) / length(mat)
})

# Write counts to output --------------------------------------------------


for (jmark in jmarks){
  print(jmark)
  mat <- mats.binfilt.lst[[jmark]]
  count.dat <- list()
  count.dat$counts <- mat
  save(count.dat, file = file.path(outdir, paste0("B6_", jmark, "_pcutoff_", pcutoff, "_binfilt_cellfilt.", Sys.Date(), ".RData")))
}


