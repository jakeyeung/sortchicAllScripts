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

# ClipLast <- function(x, jsep = "-"){
#   # B6-13W1-BM-H3K4me3-1_269 -> B6-13W1-BM-H3K4me3
#   jsplit <- strsplit(x, jsep)[[1]]
#   # remove last one
#   N <- length(jsplit) - 1
#   return(paste(jsplit[1:N], collapse = jsep))
# }
# 
# KeepLast <- function(x, jsep = "-"){
#   # B6-13W1-BM-H3K4me3-1_269 -> 1_269
#   jsplit <- strsplit(x, jsep)[[1]]
#   # remove last one
#   N <- length(jsplit)
#   return(jsplit[[N]])
# }


# Constants ---------------------------------------------------------------

# same as paper 
pcellsize <- 0.95
pfrac <- 0.95
pcutoff.bin <- 0.95


# Get quality control  ----------------------------------------------------

jmark <- "H3K4me3"

jpat <- paste0("*.", jmark, ".*.csv")
jpatall <- paste0("*.csv")
jpat.mats <- paste0("*.", jmark, ".*.rds")
jpatall.mats <- paste0("*.rds")
# jpat.mats <- paste0("*.", jmark, ".*.rds")

infs.rz <- list.files(path="/Users/yeung/data/scchic/from_cluster/count_mats_ZellerRawDataZF_all/RZcounts", pattern = jpat, full.names = TRUE)
# infs.rz <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/LHcounts", pattern = jpat, full.names = TRUE)
infs.mats <- list.files(path="/Users/yeung/data/scchic/from_cluster/count_mats_ZellerRawDataZF_all/countTables", pattern = jpat.mats, full.names = TRUE)
# infs.mats <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/countTables", pattern = jpat.mats, full.names = TRUE)
infs.rz.all <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_ZellerRawDataZF_all/RZcounts", pattern = jpatall, full.names = TRUE)
infs.mats.all <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_ZellerRawDataZF_all/countTables", pattern = jpatall.mats, full.names = TRUE)

# name by mark 
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
infs.mats.all.name <- sapply(infs.mats.all, function(x) strsplit(x, "-")[[1]][[4]], USE.NAMES = FALSE)
names(infs.mats.all) <- infs.mats.all.name

assertthat::assert_that(length(infs.rz) == length(infs.mats))

# Plot counts -------------------------------------------------------------

empty.wells <- GetEmptyWells()
dats <- ReadLH.SummarizeTA(infs.rz.all, remove.nones = FALSE) %>%
  rowwise() %>% 
  mutate(TA.count = ifelse(is.na(TA.count), 0, TA.count), 
         TA.frac = ifelse(is.na(TA.frac), 0, TA.frac),
         experi = ClipLast(samp), 
         plate = strsplit(KeepLast(samp), "_")[[1]][[1]],
         cellindx = paste("cell", strsplit(KeepLast(samp), "_")[[1]][[2]], sep = ""),
         is.empty = cellindx %in% empty.wells,
         mark = strsplit(experi, "-"))

dats <- dats %>%
  ungroup() %>%
  mutate(good.cellsize = ifelse(total.count > quantile(total.count, probs = (1-pcellsize)), TRUE, FALSE),
         good.frac = ifelse(TA.frac > quantile(TA.frac, probs = (1-pfrac)), TRUE, FALSE),
         good.cell = as.logical(good.cellsize * good.frac))

m.qc <- ggplot(dats, aes(x = log10(total.count), y = TA.frac, shape = is.empty, size = is.empty, color = good.cell)) + geom_point() + 
  facet_wrap(experi ~ plate) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load mats ---------------------------------------------------------------

mats <- lapply(infs.mats.all, function(inf){
  return(readRDS(inf))
})

rows.common <- Reduce(intersect, lapply(mats, function(mat) rownames(mat)))

# keep only from standard chromos
chromos.all <- unique(sapply(rows.common, function(x) strsplit(x, ":")[[1]][[1]]))
chromos.keep <- paste("chr", seq(25), sep = "")

bins.chromofilt.i <- which(!startsWith(rows.common, "chrK"))

rows.common.filt <- rows.common[bins.chromofilt.i]

cells.keep <- subset(dats, good.cell)$samp


# Merge mats --------------------------------------------------------------

mats.filt.merge <- lapply(jmarks, function(jmark){
  mats.sublst.i <- which(names(mats) == jmark)
  mats.sublst <- mats[mats.sublst.i]
  mats.sublst.filt <- lapply(mats.sublst, function(x){
    rows.i <- rownames(x) %in% rows.common.filt
    cols.i <- colnames(x) %in% cells.keep
    return(x[rows.i, cols.i])
  })
  return(do.call(cbind, mats.sublst.filt))
}) 

# mats.filt.merge <- do.call(cbind, mats.filt.merge)


# Find correlated bins  ---------------------------------------------------

# grep for each mats and then get bins for each mark
bin.means.long <- lapply(jmarks, function(jmark){
  mats.tmp <- mats.filt.merge[grep(jmark, names(mats.filt.merge))]
  x <- Matrix::rowMeans(do.call(cbind, mats.tmp))
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
  return(dat)
}) %>%
  bind_rows()

bin.wide.nofilt <- spread(bin.means.long, key = mark, value = bin.mean)

GGally::ggpairs(log10(bin.wide.nofilt %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic()

bin.sum <- bin.means.long %>%
  group_by(mark) %>%
  mutate(bin.rank = rank(bin.mean) / length(bin.mean)) %>%
  group_by(coord) %>%
  summarise(bin.rank.mean = mean(bin.rank),
            bin.rank.min = min(bin.rank),
            bin.rank.max = max(bin.rank))

# remove bad bins?
bins.high <- subset(bin.sum, bin.rank.mean >= pcutoff.bin)$coord
bins.low <- subset(bin.sum, bin.rank.mean <= (1-pcutoff.bin))$coord


# write PDF
outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-ZF-All"
dir.create(outdir)

pdf(file.path(outdir, paste0("qc_ZF_All.AllMarks", Sys.Date(), ".pdf")), useDingbats = FALSE)
print(m.qc)
GGally::ggpairs(log10(bin.wide.nofilt %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic()
dev.off()

# do it per mark
print(unique(dats$experi))
jdate <- "2019-11-23"
for (jmark in jmarks){
  print(jmark)
  mats.filt.merge.mark <- mats.filt.merge[[jmark]]
  
  # save to output 
  outf.all <- file.path(outdir, paste0("PZ-ZF-All_Merged.", jmark, ".", jdate, ".rds"))
  outf.unenriched <- file.path(outdir, paste0("PZ-ZF-All_Unenriched.", jmark, ".", jdate, ".rds"))
  outf.stemcells <- file.path(outdir, paste0("PZ-ZF-All_Stemcells.", jmark, ".", jdate, ".rds"))
  outf.unenrichedxstemcells <- file.path(outdir, paste0("PZ-ZF-All_UnenrichedXStemCells.", jmark, ".", jdate, ".rds"))
  
  # unenriched only
  jgrep.unenriched <- paste0("PZ-ChIC-ZFWKM-", jmark)
  # jgrep.stemcells <- paste0("PZ-ChIC-ZFWKM-", jmark)
  jgrep.stemcells <- paste0("PZ-ChIC-ZFWKMCD41plus-", jmark)
  
  cols.unenriched <- sapply(colnames(mats.filt.merge.mark), startsWith, prefix = jgrep.unenriched)
  mats.filt.unenriched <- mats.filt.merge.mark[, cols.unenriched]
  print(dim(mats.filt.unenriched))
  # assertthat::assert_that(ncol(mats.filt.unenriched) > 0)
  
  # stem cell only 
  cols.stemcells <- sapply(colnames(mats.filt.merge.mark), startsWith, prefix = jgrep.stemcells)
  mats.filt.stemcells <- mats.filt.merge.mark[, cols.stemcells]
  print(dim(mats.filt.stemcells))
  # assertthat::assert_that(ncol(mats.filt.stemcells) > 0)
  
  # # unenriched and stem cell
  # cols.unenrichedxstemcells <- grepl(paste0("^", jgrep.unenriched), colnames(mats.filt.merge.mark)) | grepl(paste0("^", jgrep.stemcells), colnames(mats.filt.merge.mark))
  # mats.filt.unenrichedxstemcells <- mats.filt.merge.mark[, cols.unenrichedxstemcells]
  # print(dim(mats.filt.unenrichedxstemcells))
  # # assertthat::assert_that(ncol(mats.filt.unenrichedxstemcells) > 0)
  # # assertthat::assert_that(ncol(mats.filt.unenrichedxstemcells) == (ncol(mats.filt.unenriched) + ncol(mats.filt.stemcells)))
  if (ncol(mats.filt.unenriched) > 0 & ncol(mats.filt.stemcells) > 0){
    mats.filt.unenrichedxstemcells <- cbind(mats.filt.unenriched, mats.filt.stemcells)
  } else {
    mats.filt.unenrichedxstemcells <- Matrix(nrow = 0, ncol = 0)
  }
  
  # write to outputs
  if (!file.exists(outf.all)){
    if (ncol(mats.filt.merge.mark) > 0){
      print(paste("Writing to", outf.all))
      saveRDS(mats.filt.merge.mark, file = outf.all)
    } else {
      print(paste("Skipping", outf.all, "empty, skipping"))
    }
  }
  if (!file.exists(outf.unenriched)){
    if (ncol(mats.filt.unenriched) > 0){
      print(paste("Writing to", outf.unenriched))
      saveRDS(mats.filt.unenriched, file = outf.unenriched)
    } else {
      print(paste("Skipping", outf.unenriched, "empty, skipping"))
    }
  }
  if (!file.exists(outf.stemcells)){
    if (ncol(mats.filt.stemcells) > 0){
      print(paste("Writing to", outf.stemcells))
      saveRDS(mats.filt.stemcells, file = outf.stemcells)
    } else {
      print(paste("Skipping", outf.stemcells, "empty, skipping"))
    }
  }
  if (!file.exists(outf.unenrichedxstemcells)){
    if (ncol(mats.filt.unenrichedxstemcells) > 0){
      print(paste("Writing to", outf.unenrichedxstemcells))
      saveRDS(mats.filt.unenrichedxstemcells, file = outf.unenrichedxstemcells)
    } else {
      print(paste("Skipping", outf.unenrichedxstemcells, "empty, skipping"))
    }
  }
}






