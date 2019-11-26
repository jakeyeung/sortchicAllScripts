# Jake Yeung
# Date of Creation: 2019-11-21
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/1-explore_unenriched_linneg_stemcell_all_merged.R
# Rerun on same pipeline to see if cells "mix" better

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(scchicFuncs)

# Get enriched -------------------------------------------------------------

# read count mats per mark
jmark <- "H3K9me3"
# bad.experi <- "B6-13W1-BM-H3K4me3"
bad.experi <- ""

outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_Bl6-BM-Redo"

# look at unenriched
indir <- "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg/ZellerRawDataB6/countTables"
indir.rz <- "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg/ZellerRawDataB6/RZcounts"

infs <- list.files(indir, pattern = paste0("*.", jmark, ".*.countTable.rds"), full.names = TRUE)
infs.rz <- list.files(indir.rz, pattern = paste0("*.", jmark, ".*.csv"), full.names = TRUE)


dats <- lapply(infs.rz, function(inf){
  dat <- ReadRZ(inf)
  # make long
  dat.TA <- tidyr::gather(dat %>% filter(V1 == "TA_start"), key = "samp", value = "count", -V1) %>%
    dplyr::rename(TA.count = count) %>%
    dplyr::select(-V1)
  dat.Total <- tidyr::gather(dat %>% filter(V1 == "total"), key = "samp", value = "count", -V1) %>%
    dplyr::rename( total.count = count) %>%
    dplyr::select(-V1)
  # merge the two
  dat.long <- left_join(dat.TA, dat.Total, by = c("samp"))
}) %>%
  bind_rows() %>%
  mutate(TA.frac = TA.count / total.count)

  
# Get Index and Experi
dats$experi <- sapply(dats$samp, function(x) strsplit(x, "_")[[1]][[1]])
dats$indx <- sapply(dats$experi, function(x){
  jsplit <- strsplit(x, "-")[[1]]
  return(jsplit[length(jsplit) - 3])
})
dats$is.stem.enriched <- sapply(dats$experi, function(x){
  return(grepl("stem-cell", x))
})
dats$cell <- sapply(dats$samp, function(x){
  return(paste0("cell", strsplit(x, "_")[[1]][[2]]))
})

# check samp names are unique
assertthat::assert_that(length(which(duplicated(subset(dats)$samp))) == 0)

ggplot(dats, aes(x = log10(total.count), y = TA.frac)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~experi)

empty.wells <- GetEmptyWells()

dats$is.empty <- sapply(dats$cell, function(x) x %in% empty.wells)


# set a cutoff
jcutoff <- 2.5
ta.cutoff <- 0.6

dats <- dats %>%
  rowwise() %>%
  mutate(good.cell = (TA.frac >= ta.cutoff & log10(total.count) >= jcutoff))

ggplot(dats, aes(x = log10(total.count), y = TA.frac, color = good.cell, shape = is.empty)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~experi) +
  geom_vline(xintercept = jcutoff) +
  geom_hline(yintercept = ta.cutoff)


good.cells <- subset(dats, experi != bad.experi & good.cell)$cell




# Get mats ----------------------------------------------------------------

names(infs) <- sapply(infs, function(x){
  strsplit(basename(x), "\\.")[[1]][[1]]
})

prefix <- unique(sapply(names(infs), function(x) gsub("-[0-9]$", "", x)))

mats <- lapply(infs, function(inf){
  return(readRDS(inf))
}) 

# remove bad plates
mats.filt <- lapply(mats, function(x){
  if (nrow(x) < 50000){
    return(NULL)
  } else {
    return(x)
  }
}) %>%
  purrr::compact()

# filter common rows for now
common.rows <- lapply(mats.filt, function(x){
  return(rownames(x))
}) %>%
  Reduce(intersect, .)

mats.filt <- lapply(mats.filt, function(x){
  rows.i <- rownames(x) %in% common.rows
  x[rows.i, ]
})

mats.merge <- do.call(cbind, mats.filt)

# Filter good cells -------------------------------------------------------

outf <- file.path(outdir, paste0(prefix, ".good_cells_filt.rds"))

# mats.merge.filt <- 










