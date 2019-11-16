# Jake Yeung
# Date of Creation: 2019-11-13
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/H3K27me3_QC_get_good_cells_bins.R
# Get good cells and bins quality control K27me3

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)

library(Matrix)

jcutoff <- 2.5
ta.cutoff <- 0.5

# Load just H3K27me3 ------------------------------------------------------

# shoudl be 4 plates 
dir.rz <- "/Users/yeung/data/scchic/from_cluster/oud3985/RZcounts"
dir.count <- "/Users/yeung/data/scchic/from_cluster/oud3985/countTables"

gstr <- "*ZFWKM-H3K27me3.*.csv"
gstr2 <- "*ZFWKM-H3K27me3.*.rds"

(infs.rz <- list.files(dir.rz, pattern = gstr, all.files = TRUE, full.names = TRUE))
(infs.count <- list.files(dir.count, pattern = gstr2, all.files = TRUE, full.names = TRUE))

bname <- unique(sapply(infs.count, function(x) gsub(pattern = "-[0-9]$",replacement = "", x = strsplit(basename(x), split = "\\.")[[1]][[1]])))

# Read tables -------------------------------------------------------------

dat.rz <- lapply(infs.rz, function(inf){
  dat <- ReadRZ(inf, remove.nones = TRUE)
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
  rowwise() %>%
  mutate(TA.count = ifelse(is.na(TA.count), 0, TA.count), 
         TA.frac = TA.count / total.count,
         experi = strsplit(samp, "_")[[1]][[1]])

m.qc <- ggplot(dat.rz, aes(x = log10(total.count), y = TA.frac)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = ta.cutoff) + geom_vline(xintercept = jcutoff) + facet_wrap(~experi)
print(m.qc)

dat.rz <- dat.rz %>%
  rowwise() %>%
  mutate(good.cell = log10(total.count) > jcutoff & TA.frac > ta.cutoff)

# show empty wells
empty.wells <- GetEmptyWells()

dat.rz <- dat.rz %>%
  rowwise() %>%
  mutate(cell = paste("cell", strsplit(samp, split = "_")[[1]][[2]], sep = ""),
         is.empty = cell %in% empty.wells)

m.qc <- ggplot(dat.rz, aes(x = log10(total.count), y = TA.frac, color = good.cell, shape = is.empty, size = is.empty)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = ta.cutoff) + geom_vline(xintercept = jcutoff) + facet_wrap(~experi)
print(m.qc)

cells.keep <- subset(dat.rz, good.cell)$samp

# Get good bins -----------------------------------------------------------

# take from existing LDA 
jmark <- "H3K4me3"
jbin <- "FALSE"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

terms.keep <- rownames(count.mat)

# Make count tables -------------------------------------------------------

# load count tables


# keep all bins and filter them out later 
mats <- lapply(infs.count, function(inf){
  mat.tmp <- readRDS(inf)
})

lapply(mats, function(x) nnzero(x) / length(x))
lapply(mats, dim)

rnames.lst <- lapply(mats, function(x) rownames(x))

terms.keep.common <- Reduce(f = intersect, x = rnames.lst)
terms.keep.common <- terms.keep.common[terms.keep.common %in% terms.keep]

# filter bins and good cells
mats.filt <- lapply(mats, function(mat){
  rows.i <- rownames(mat) %in% terms.keep.common
  cols.i <- colnames(mat) %in% cells.keep
  return(mat[rows.i, cols.i])
})
lapply(mats.filt, dim)

mat.merge <- do.call(cbind, mats.filt)


# Write good cells  -------------------------------------------------------




# Write to output ---------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow.H3K27me3"
outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow"
dir.create(outdir)

good.cells.path <- file.path(outdir, paste0("ZF-H3K27me3.", Sys.Date(), ".csv"))
sink(file = good.cells.path)
for (jcell in cells.keep){
  cat(paste0(jcell, "\n"))
}
sink()

# saveRDS(mat.merge, file = file.path(outdir, paste0(bname, ".", Sys.Date(), ".rds")))








