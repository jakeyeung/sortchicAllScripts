# Jake Yeung
# Date of Creation: 2019-06-17
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/2-quality_control_prepare_matrix.R
# Prepare count matrix

library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(tidyr)
library(GGally)

source("scripts/Rfunctions/QCFunctions.R")
source("scripts/Rfunctions/Aux.R")

# Constants  -------------------------------------------------------------------

TAcutoff <- 0.7
countcutoff <- 3

prefix <- "_PZ-Bl6-BM-Linneg"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


bamlist.outdir <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_peak_analysis_build95", prefix)
outdir <- paste0("~/data/scchic/count_mat_binfilt_cellfilt_for_LDA", prefix)
dir.create(outdir)
dir.create(bamlist.outdir)

indir <- paste0("/Users/yeung/data/scchic/from_cluster/count_mat_from_bam", prefix)
assertthat::assert_that(dir.exists(indir))


# Load --------------------------------------------------------------------



infs.h3k4me1 <- list.files(indir, pattern = ".*H3K4me1.*.rds", full.names = TRUE)
infs.h3k4me3 <- list.files(indir, pattern = ".*H3K4me3.*.rds", full.names = TRUE)
infs.h3k27me3 <- list.files(indir, pattern = ".*H3K27me3.*.rds", full.names = TRUE)
infs.h3k9me3 <- list.files(indir, pattern = ".*H3K9me3.*.rds", full.names = TRUE)
infs.all <- list.files(indir, pattern = "*.rds", full.names = TRUE)
infs.bymark <- list(infs.h3k4me1, infs.h3k4me3, infs.h3k27me3, infs.h3k9me3)
names(infs.bymark) <- jmarks

inf.cells <- paste0("~/data/scchic/quality_control", prefix, "/good_cells.TAcutoff.", TAcutoff, ".countcutoff.", countcutoff, ".txt")
assertthat::assert_that(file.exists(inf.cells))
cells <- fread(inf.cells, header = FALSE, col.names = c("cell"))
cells.keep <- cells$cell

bin.means.h3k4me1 <- Matrix::rowMeans(LoadMats(infs.h3k4me1, cells.keep))
bin.means.h3k4me3 <- Matrix::rowMeans(LoadMats(infs.h3k4me3, cells.keep))
bin.means.h3k27me3 <- Matrix::rowMeans(LoadMats(infs.h3k27me3, cells.keep))
bin.means.h3k9me3 <- Matrix::rowMeans(LoadMats(infs.h3k9me3, cells.keep))


bin.lst <- list(bin.means.h3k4me1, bin.means.h3k4me3, bin.means.h3k27me3, bin.means.h3k9me3)
names(bin.lst) <- jmarks

rnames.all <- lapply(bin.lst, function(x) names(x))

rnames.all.common <- Reduce(intersect, rnames.all)



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



# load bad bins from previous analysis ------------------------------------

bins.correlated.dat <- fread("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_correlated_bins.2019-05-11.bed")
bins.correlated <- paste(bins.correlated.dat$V1, paste(bins.correlated.dat$V2, bins.correlated.dat$V3, sep = "="), sep = ":")

# check any systematic biases on the offset ------------------------------

bsize <- 100000
bin.sum$bin <- as.numeric(sapply(bin.sum$coord, GetStart))
bin.sum$offset <- sapply(bin.sum$bin, function(x) x %% bsize)

ggplot(bin.sum, aes(x = bin.rank.mean, group = offset)) + geom_histogram() + facet_wrap(~offset, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rnames.blfilt.corrfilt <- (bin.wide %>% dplyr::filter(!coord %in% bins.correlated))$coord


# Write correlated bins to output -----------------------------------------


for (jmark in jmarks){
  cells.keep.bymark <- data.frame(cell = paste0(grep(pattern = jmark, cells.keep, value = TRUE), ".sorted.bam"))
  fwrite(cells.keep.bymark, file = file.path(bamlist.outdir, paste0("JY_", jmark, "_bamnames.out")), col.names = FALSE, sep = "\t")
}


# Do final bin filter on matrix -------------------------------------------

mats.lst <- lapply(infs.bymark, LoadMats, cells.keep)

lapply(mats.lst, dim)
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

prefix2 <- substring(prefix, first = 2, last = nchar(prefix))
for (jmark in jmarks){
  print(jmark)
  mat <- mats.binfilt.lst[[jmark]]
  count.dat <- list()
  count.dat$counts <- mat
  save(count.dat, file = file.path(outdir, paste0(prefix2, "_", jmark, "_binfilt_cellfilt.", Sys.Date(), ".RData")))
}


# Load the previous analysis and integrate with the count dat there -------

indir.orig <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA"
infs.orig <- lapply(jmarks, function(jmark) file.path(indir.orig, paste0("B6_", jmark, "_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData")))

mats.binlist.lst.orig <- lapply(jmarks, function(jmark){
  load(infs.orig[[jmark]], v=T)
  return(count.dat$counts)
})

# concatenate with original matrix and run 
for (jmark in jmarks){
  print(jmark)
  mat.orig <- mats.binlist.lst.orig[[jmark]]
  rownames(mat.orig) <- paste("chr", rownames(mat.orig), sep="")
  mat.linneg <- mats.binfilt.lst[[jmark]]
  rows.common <- intersect(rownames(mat.orig), rownames(mat.linneg))
  mat.merged <- Matrix::cbind2(mat.orig[rows.common, ], mat.linneg[rows.common, ])
  count.dat <- list()
  count.dat$counts <- mat.merged
  save(count.dat, file = file.path(outdir, paste0(prefix2, "_", jmark, "_binfilt_cellfilt.", Sys.Date(), ".merge_with_B6.RData")))
}

