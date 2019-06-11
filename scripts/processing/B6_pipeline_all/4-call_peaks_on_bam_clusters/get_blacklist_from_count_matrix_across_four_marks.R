# Jake Yeung
# Date of Creation: 2019-04-21
# File: ~/projects/scchic/scripts/scripts_analysis/tf_activity_debugging/explore_count_matrix_filter_correlated_peaks.R
# Remove bad peaks

rm(list=ls())

tstart <- Sys.time()

library(dplyr)
library(ggplot2)
library(Matrix)
library(GGally)
library(tidyr)

source("scripts/Rfunctions/Aux.R")

args <- commandArgs(trailingOnly=TRUE)

inmain <- args[[1]]
refmark <- args[[2]]
jcutoff <- as.numeric(args[[3]])
outmain <- args[[4]]

assertthat::assert_that(is.numeric(jcutoff))
assertthat::assert_that(dir.exists(inmain))
assertthat::assert_that(dir.exists(outmain))


# Constnats ---------------------------------------------------------------

# jcutoff <- 0.99

# Load data ---------------------------------------------------------------

# refmark <- "H3K4me1"
# refmark <- "H3K4me3"
# refmark <- "H3K27me3"
marks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(marks) <- marks

# inmain <- paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.1000_build95.withchr.cells_from_bin_analysis")
dname <- paste0("peak_from_", refmark)
indir <- file.path(inmain, dname)
print(paste("Ref mark", refmark))
peak.mean.dat <- lapply(marks, function(mark){
  print(mark)
  # inf <- paste0("/Users/yeung/data/scchic/from_cluster/count_mat_explore.peak_ref_", refmark, "/PZ-BM-", mark, ".merged.NoCountThres.hiddenDomains.Robj")
  inf <- file.path(indir, paste0("PZ-BM-", mark, ".merged.NoCountThres.hiddenDomains.Robj"))
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  count.filt <- count.dat$counts
  # count.filt[!is.na(rownames(count.filt))] <- count.filt
  count.filt <- sweep(count.filt, MARGIN = 2, STATS = colSums(count.filt), FUN = "/")
  peakmean <- rowMeans(count.filt)
  dat.out <- data.frame(peak = rownames(count.filt), count.mean = peakmean, mark = mark, stringsAsFactors = FALSE)
  return(dat.out) 
}) %>% bind_rows()

peak.mean.dat$start <- as.numeric(sapply(peak.mean.dat$peak, GetStart))
peak.mean.dat$end <- as.numeric(sapply(peak.mean.dat$peak, GetEnd))
peak.mean.dat$bsize <- peak.mean.dat$end - peak.mean.dat$start

peak.mean.dat <- peak.mean.dat %>%
  group_by(mark) %>%
  mutate(count.mean.norm = count.mean / bsize,
         p = rank(count.mean.norm) / length(count.mean.norm))

peak.mean.sum <- peak.mean.dat %>%
  group_by(peak, bsize) %>%
  summarise(pmean = mean(p)) %>%
  arrange(desc(pmean))

# check cutoffs
subset(peak.mean.sum, pmean < jcutoff)

# keep 

peaks.high <- as.character((peak.mean.sum %>% filter(pmean > jcutoff))$peak)
peaks.keep <- as.character((peak.mean.sum %>% filter(pmean < jcutoff))$peak)
peaks.lowfilt <- as.character((peak.mean.sum %>% filter(pmean > 0.4))$peak)
peaks.plot <- intersect(peaks.keep, peaks.lowfilt)
print(length(peaks.plot))

peak.mat <- spread(data = peak.mean.dat %>% dplyr::select(-p), key = mark, value = count.mean)

# write a new blacklist
bl.dat <- data.frame(chromo = sapply(peaks.high, GetChromo),
                     start = sapply(peaks.high, GetStart),
                     end = sapply(peaks.high, GetEnd), 
                     stringsAsFactors = FALSE) %>% arrange(chromo, start, end)

# mm10_BM_H3K4me1_blacklist_QuantileCutoff-0.99.bed.gz

outfname <- paste0("mm10_BM_", refmark, "_blacklist_QuantileCutoff-", jcutoff, ".bed")
outf <- file.path(outmain, outfname)
data.table::fwrite(bl.dat, file = outf, sep = "\t", col.names = FALSE)
system(paste0("gzip ", outf))

print(tstart - Sys.time())
