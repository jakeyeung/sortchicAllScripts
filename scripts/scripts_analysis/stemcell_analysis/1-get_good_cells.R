# Jake Yeung
# Date of Creation: 2019-09-29
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/1-get_good_cells.R
# Get good cells

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)


# Functions ---------------------------------------------------------------

ReadRZ <- function(inf){
  dat.tmp <- as.data.frame(fread(inf, sep = ",", header = TRUE))
  colnames(dat.tmp)[[1]] <- dat.tmp[1, 1]
  dat.tmp <- dat.tmp[-1, ]
  dat.mat <- as.matrix(dat.tmp[, -1])
  dat.mat.sum <- colSums(dat.mat, na.rm = TRUE)
  dat.mat.ta <- subset(dat.tmp, recognizedSequence == "TA", select = -recognizedSequence)  # I think Buys calls it AT for some reason here? Need to check actually, might be an old bam file
  dat.merged <- rbind(dat.mat.ta, dat.mat.sum)
  # add first column
  dat.merged.annot <- cbind(data.frame(V1 = c("TA_start", "total")), dat.merged)  # call it V1 and rownames same as Buys TA_obs_per_cell.csv
  return(dat.merged.annot)
}


# Outputs -----------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells"

# Define files ------------------------------------------------------------



# first the RZ
indir.rz <- "/Users/yeung/data/scchic/from_cluster/oud3700/RZcounts"
infs.rz <- list.files(indir.rz, pattern = "*RZ_counts.csv", full.names = TRUE)

# count tables
indir.count <- "/Users/yeung/data/scchic/from_cluster/oud3700/countTables"
infs.count <- list.files(indir.count, pattern = "*countTable.rds", full.names = TRUE)



# Load data ---------------------------------------------------------------


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
 

# Load good H3K4me1 bins and write to output  -----------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1"); names(jmarks) <- jmarks
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
tm.result.lst <- lapply(infs, LoadGetTmResult)

jterms <- colnames(tm.result.lst$H3K4me1$terms)

infs.count.bnames <- sapply(infs.count, basename)

mats <- lapply(infs.count, function(inf){
  mat.tmp <- readRDS(inf)
  rows.keep <- which(rownames(mat.tmp) %in% jterms)
  return(mat.tmp[rows.keep, ])
})

lapply(mats, function(x) nnzero(x) / length(x))
lapply(mats, dim)

# get common rows
rnames <- lapply(mats, rownames)
rnames.common <- Reduce(intersect, rnames)

# do one more filter and write to output
mats.filt <- lapply(mats, function(x){
  rows.keep <- which(rownames(x) %in% rnames.common)
  return(x[rows.keep, ])
})
lapply(mats.filt, dim)

names(mats.filt) <- infs.count.bnames


# Write to output ---------------------------------------------------------

# do one with all
mats.merged.all <- do.call(cbind, mats.filt)
print(dim(mats.merged.all))

# do one with only nonenriched
mats.keep <- !grepl("stem-cells", infs.count.bnames)
mats.merged.nonenriched <- do.call(cbind, mats.filt[mats.keep])
print(dim(mats.merged.nonenriched))

# do one with only enriched
mats.keep.enrich <- grepl("stem-cells", infs.count.bnames)
mats.merged.enriched <- do.call(cbind, mats.filt[mats.keep.enrich])
print(dim(mats.merged.enriched))

# do one with only enriched, removing Index3 (all cell cycles)
mats.keep.enrich.names <- grep("stem-cells", infs.count.bnames, value = TRUE)
mats.keep.enrich.noCellCycle <- grepl("Index3", mats.keep.enrich.names)
mats.merged.enriched.noCellCycle <- do.call(cbind, mats.filt[mats.keep.enrich.noCellCycle])
print(dim(mats.merged.enriched.noCellCycle))


mats.lst <- c(mats.merged.all, mats.merged.nonenriched, mats.merged.enriched, mats.merged.enriched.noCellCycle)
names(mats.lst) <- c("matsMergedAll", "matsMergedNonenriched", "matsMergedEnriched", "matsMergedEnrichedNoCellCycle")

# Write output names ------------------------------------------------------

jmark <- "H3K4me1"
lapply(names(mats.lst), function(jname){
  print(jname)
  mats.merge <- mats.lst[[jname]]
  print(dim(mats.merge))
  count.dat <- list()
  count.dat$counts <- mats.merge
  save(count.dat, file = file.path(outdir, paste0("PZ-Bl6-BM-StemCells_", jmark, "_", jname, "_", Sys.Date(), ".RData")))
})



