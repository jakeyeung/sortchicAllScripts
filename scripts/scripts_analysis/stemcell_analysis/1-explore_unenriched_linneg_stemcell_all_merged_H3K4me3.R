# Jake Yeung
# Date of Creation: 2019-11-22
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/5-analyze_stem_cells_merged.R
# Analyze everything together

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)
library(topicmodels)
library(Matrix)

ClipLast <- function(x, jsep = "-"){
  # B6-13W1-BM-H3K4me3-1_269 -> B6-13W1-BM-H3K4me3
  jsplit <- strsplit(x, jsep)[[1]]
  # remove last one
  N <- length(jsplit) - 1
  return(paste(jsplit[1:N], collapse = jsep))
}

KeepLast <- function(x, jsep = "-"){
  # B6-13W1-BM-H3K4me3-1_269 -> 1_269
  jsplit <- strsplit(x, jsep)[[1]]
  # remove last one
  N <- length(jsplit)
  return(jsplit[[N]])
}


# Constants ---------------------------------------------------------------

# same as paper 
pcellsize <- 0.85
pfrac <- 0.95


# Get quality control  ----------------------------------------------------

jmark <- "H3K4me3"

jpat <- paste0("*.", jmark, ".*.csv")
jpat.mats <- paste0("*.", jmark, ".*.rds")
infs.rz <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/LHcounts", pattern = jpat, full.names = TRUE)
infs.mats <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/countTables", pattern = jpat.mats, full.names = TRUE)

assertthat::assert_that(length(infs.rz) == length(infs.mats))

# Plot counts -------------------------------------------------------------

empty.wells <- GetEmptyWells()
dats <- ReadLH.SummarizeTA(infs.rz, remove.nones = FALSE) %>%
  rowwise() %>% 
  mutate(TA.count = ifelse(is.na(TA.count), 0, TA.count), 
         TA.frac = ifelse(is.na(TA.frac), 0, TA.frac),
         experi = ClipLast(samp), 
         plate = strsplit(KeepLast(samp), "_")[[1]][[1]],
         cellindx = paste("cell", strsplit(KeepLast(samp), "_")[[1]][[2]], sep = ""),
         is.empty = cellindx %in% empty.wells)

dats <- dats %>%
  ungroup() %>%
  mutate(good.cellsize = ifelse(total.count > quantile(total.count, probs = (1-pcellsize)), TRUE, FALSE),
         good.frac = ifelse(TA.frac > quantile(TA.frac, probs = (1-pfrac)), TRUE, FALSE),
         good.cell = as.logical(good.cellsize * good.frac))

m.qc <- ggplot(dats, aes(x = log10(total.count), y = TA.frac, shape = is.empty, size = is.empty, color = good.cell)) + geom_point() + 
  facet_wrap(experi ~ plate) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Get good bins from previous analysis ------------------------------------

inf.orig <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.", jmark, ".stringent_filter.RData")
assertthat::assert_that(file.exists(inf.orig))
load(inf.orig, v=T)

bins.keep <- out.objs$out.lda@terms
cells.keep <- subset(dats, good.cell)$samp

# Load mats ---------------------------------------------------------------

mats <- lapply(infs.mats, function(inf){
  return(readRDS(inf))
})

rows.common <- Reduce(intersect, lapply(mats, function(mat) rownames(mat)))

rows.common.filt <- rows.common[rows.common %in% bins.keep]



# Merge mats --------------------------------------------------------------

mats.filt.merge <- lapply(mats, function(x){
  rows.i <- rownames(x) %in% rows.common.filt
  cols.i <- colnames(x) %in% cells.keep
  x[rows.i, cols.i]
}) 

mats.filt.merge <- do.call(cbind, mats.filt.merge)

# write PDF
outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All"

pdf(file.path(outdir, paste0("qc_BM_All.", jmark, ".", Sys.Date(), ".pdf")), useDingbats = FALSE)
print(m.qc)
dev.off()

# save to output 
jdate <- "2019-11-22"
outf.all <- file.path(outdir, paste0("PZ-Bl6-BM-All_Merged.", jmark, ".", jdate, ".rds"))
outf.unenriched <- file.path(outdir, paste0("PZ-Bl6-BM-All_Unenriched.", jmark, ".", jdate, ".rds"))
outf.linneg <- file.path(outdir, paste0("PZ-Bl6-BM-All_Linneg.", jmark, ".", jdate, ".rds"))
outf.stemcells <- file.path(outdir, paste0("PZ-Bl6-BM-All_Stemcells.", jmark, ".", jdate, ".rds"))
outf.unenrichedxlinneg <- file.path(outdir, paste0("PZ-Bl6-BM-All_UnenrichedXLinNeg.", jmark, ".", jdate, ".rds"))
outf.unenrichedxstemcells <- file.path(outdir, paste0("PZ-Bl6-BM-All_UnenrichedXStemCells.", jmark, ".", jdate, ".rds"))

print(unique(dats$experi))

# unenriched only
cols.unenriched <- sapply(colnames(mats.filt.merge), startsWith, prefix = "B6-13W1-BM-H3K4me3")
mats.filt.unenriched <- mats.filt.merge[, cols.unenriched]
print(dim(mats.filt.unenriched))
assertthat::assert_that(ncol(mats.filt.unenriched) > 0)

# linneg only
cols.linneg <- sapply(colnames(mats.filt.merge), startsWith, prefix = "PZ-Bl6-BM-Linneg-H3K4me3")
mats.filt.linneg <- mats.filt.merge[, cols.linneg]
print(dim(mats.filt.linneg))
assertthat::assert_that(ncol(mats.filt.linneg) > 0)

# stem cell only 
cols.stemcells <- sapply(colnames(mats.filt.merge), startsWith, prefix = "PZ-ChIC-B6BMSC-")
mats.filt.stemcells <- mats.filt.merge[, cols.stemcells]
print(dim(mats.filt.stemcells))
assertthat::assert_that(ncol(mats.filt.stemcells) > 0)

# unenriched and linneg
# cols.unenriched <- sapply(colnames(mats.filt.merge), startsWith, prefix = "B6-13W1-BM-H3K4me3")
cols.unenrichedxlinneg <- grepl("^B6-13W1-BM-H3K4me3", colnames(mats.filt.merge)) | grepl("^B6-13W1-BM-H3K4me3", colnames(mats.filt.merge))
mats.filt.unenrichedxlinneg <- mats.filt.merge[, cols.unenrichedxlinneg]
print(dim(mats.filt.unenrichedxlinneg))
assertthat::assert_that(ncol(mats.filt.unenrichedxlinneg) > 0)


# unenriched and stem cell
cols.unenrichedxstemcells <- grepl("^B6-13W1-BM-H3K4me3", colnames(mats.filt.merge)) | grepl("^PZ-ChIC-B6BMSC-", colnames(mats.filt.merge))
mats.filt.unenrichedxstemcells <- mats.filt.merge[, cols.unenrichedxstemcells]
print(dim(mats.filt.unenrichedxstemcells))
assertthat::assert_that(ncol(mats.filt.unenrichedxstemcells) > 0)


# write to outputs
if (!file.exists(outf.all)){
  saveRDS(mats.filt.merge, file = outf.all)
}
if (!file.exists(outf.unenriched)){
  saveRDS(mats.filt.unenriched, file = outf.unenriched)
}
if (!file.exists(outf.linneg)){
  saveRDS(mats.filt.linneg, file = outf.linneg)
}
if (!file.exists(outf.stemcells)){
  saveRDS(mats.filt.stemcells, file = outf.stemcells)
}
if (!file.exists(outf.unenrichedxlinneg)){
  saveRDS(mats.filt.unenrichedxlinneg, file = outf.unenrichedxlinneg)
}
if (!file.exists(outf.unenrichedxstemcells)){
  saveRDS(mats.filt.unenrichedxstemcells, file = outf.unenrichedxstemcells)
}





