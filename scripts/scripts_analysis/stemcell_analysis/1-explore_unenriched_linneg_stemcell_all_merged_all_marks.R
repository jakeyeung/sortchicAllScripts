# Jake Yeung
# Date of Creation: 2019-11-23
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/1-explore_unenriched_linneg_stemcell_all_merged_all_marks.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)
library(topicmodels)
library(Matrix)


library(stringr)
GetMarkFromStr <- function(x){
  jpat <- "H3K[0-9]*.me[0-9]"
  jmatch <- stringr::str_match(x, jpat)
  assertthat::assert_that(nrow(jmatch) == 1 & ncol(jmatch) == 1)
  return(jmatch[[1]])
}


# Constants ---------------------------------------------------------------

# same as paper 
pcellsize <- 0.95
pfrac <- 0.95
pcutoff.bin <- 0.95


# Get quality control  ----------------------------------------------------

jpatall <- paste0("*.csv")
jpatall.mats <- paste0("*.rds")

infs.rz.all <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/LHcounts", pattern = jpatall, full.names = TRUE)
infs.mats.all <- list.files(path = "/Users/yeung/data/scchic/from_cluster/count_mats_B6BM_and_LinNeg_StemCells/ZellerRawDataB6_all/countTables", pattern = jpatall.mats, full.names = TRUE)

# name by mark 
# H3K4me3 is done separately because we use a more stringent cutoff??
# jmarks <- c("H3K4me1", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
infs.mats.all.name <- sapply(infs.mats.all, function(x) GetMarkFromStr(basename(x)))
names(infs.mats.all) <- infs.mats.all.name

assertthat::assert_that(length(infs.rz.all) == length(infs.mats.all))




# Plot counts -------------------------------------------------------------

empty.wells <- GetEmptyWells()
dats <- ReadLH.SummarizeTA(infs.rz.all, remove.nones = FALSE) %>%
  rowwise() %>% 
  mutate(TA.count = ifelse(is.na(TA.count), 0, TA.count), 
         TA.frac = ifelse(is.na(TA.frac), 0, TA.frac),
         experi = ClipLast(samp), 
         plate = strsplit(KeepLast(samp), "_")[[1]][[1]],
         experiplate = paste(experi, plate, sep = "_"), 
         cellindx = paste("cell", strsplit(KeepLast(samp), "_")[[1]][[2]], sep = ""),
         is.empty = cellindx %in% empty.wells,
         mark = GetMarkFromStr(experi))

dats <- dats %>%
  ungroup() %>%
  mutate(good.cellsize = ifelse(total.count > quantile(total.count, probs = (1-pcellsize)), TRUE, FALSE),
         good.frac = ifelse(TA.frac > quantile(TA.frac, probs = (1-pfrac)), TRUE, FALSE),
         good.cell = as.logical(good.cellsize * good.frac))

m.qc <- ggplot(dats, aes(x = log10(total.count), y = TA.frac, shape = is.empty, size = is.empty, color = good.cell)) + geom_point() + 
  facet_wrap(experi ~ plate) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# remove bad experiments
min.good.cells <- 300
good.cells.by.plate <- dats %>%
  group_by(experiplate) %>%
  summarise(n.goodcells = length(which(good.cell))) %>%
  arrange(n.goodcells)

bad.plates <- subset(good.cells.by.plate, n.goodcells < min.good.cells)

print(paste("removing bad plates:", bad.plates$experiplate))

dats <- subset(dats, !experiplate %in% bad.plates)

cells.keep <- subset(dats, good.cell)$samp  # filtered out bad plates already 

# Load mats ---------------------------------------------------------------

mats <- lapply(infs.mats.all, function(inf){
  jmat <- readRDS(inf) 
  # filter out bad cells
  cols.i <- which(colnames(jmat) %in% cells.keep)
  if (length(cols.i) > 0){
    return(jmat[, cols.i])
  } else{
    return(NULL)
  }
})

# remove NULLs in list
mats <- purrr::compact(mats)


# merge by mark
mats.merge <- lapply(jmarks, function(jmark){
  mats.sublst.i <- which(names(mats) == jmark)
  mats.filt <- mats[mats.sublst.i]
  
  # get temporary common rows (AFTER removing bad cells)
  rows.common.temp <- Reduce(intersect, lapply(mats.filt, rownames))
  mats.filt.filt <- lapply(mats.filt, function(x){
    rows.i <- which(rownames(x) %in% rows.common.temp)
    return(x[rows.i, ])
  }) %>%
    do.call(cbind, .)
})

print(lapply(mats.merge, dim))

rows.common <- Reduce(intersect, lapply(mats, function(mat) rownames(mat)))



# remake good bins
bins.keep.lst <- lapply(jmarks, function(jmark){
  inf.orig <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.", jmark, ".RData")
  assertthat::assert_that(file.exists(inf.orig))
  load(inf.orig, v=T)
  print(out.objs$out.lda@k)
  return(out.objs$out.lda@terms)
})

bins.keep <- Reduce(intersect, bins.keep.lst)

rows.common.filt <- rows.common[which(rows.common %in% bins.keep)]

print(length(rows.common.filt))

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


# write PDF
outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks"
dir.create(outdir)

pdf(file.path(outdir, paste0("qc_ZF_All.AllMarks", Sys.Date(), ".pdf")), useDingbats = FALSE)
print(m.qc)
# GGally::ggpairs(log10(bin.wide.nofilt %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic()
dev.off()

# do it per mark
print(unique(dats$experi))
jdate <- Sys.Date()

jgrep.strs <- lapply(jmarks, function(jmark){
  jgrep.unenriched <- paste0("B6-13W1-BM-", jmark)
  jgrep.linneg <- paste0("PZ-Bl6-BM-Linneg-", jmark)
  jgrep.stemcells <- ifelse(jmark == "H3K4me1", paste0("PZ-ChIC-Bl6-BM-stem-cells-", jmark), paste0("PZ-ChIC-B6BMSC-", jmark))
  return(list(unenriched = jgrep.unenriched, linneg = jgrep.linneg, stemcells = jgrep.stemcells))
})

for (jmark in jmarks){
  print(jmark)
  mats.filt.merge.mark <- mats.filt.merge[[jmark]]
  outf.all <- file.path(outdir, paste0("PZ-Bl6-BM-All_Merged.", jmark, ".", jdate, ".rds"))
  outf.unenriched <- file.path(outdir, paste0("PZ-Bl6-BM-All_Unenriched.", jmark, ".", jdate, ".rds"))
  outf.linneg <- file.path(outdir, paste0("PZ-Bl6-BM-All_Linneg.", jmark, ".", jdate, ".rds"))
  outf.stemcells <- file.path(outdir, paste0("PZ-Bl6-BM-All_Stemcells.", jmark, ".", jdate, ".rds"))
  outf.unenrichedxlinneg <- file.path(outdir, paste0("PZ-Bl6-BM-All_UnenrichedXLinNeg.", jmark, ".", jdate, ".rds"))
  outf.unenrichedxstemcells <- file.path(outdir, paste0("PZ-Bl6-BM-All_UnenrichedXStemCells.", jmark, ".", jdate, ".rds"))
  
  print(unique(dats$experi))
  
  jgrep.unenriched <- jgrep.strs[[jmark]]$unenriched
  jgrep.linneg <- jgrep.strs[[jmark]]$linneg
  jgrep.stemcells <- jgrep.strs[[jmark]]$stemcells
  
  # unenriched only
  cols.unenriched <- sapply(colnames(mats.filt.merge.mark), startsWith, prefix = jgrep.unenriched)
  mats.filt.unenriched <- mats.filt.merge.mark[, cols.unenriched]
  print(dim(mats.filt.unenriched))
  # assertthat::assert_that(ncol(mats.filt.unenriched) > 0)
  
  # linneg only
  cols.linneg <- sapply(colnames(mats.filt.merge.mark), startsWith, prefix = jgrep.linneg)
  mats.filt.linneg <- mats.filt.merge.mark[, cols.linneg]
  print(dim(mats.filt.linneg))
  # assertthat::assert_that(ncol(mats.filt.linneg) > 0)
  
  # stem cell only 
  cols.stemcells <- sapply(colnames(mats.filt.merge.mark), startsWith, prefix = jgrep.stemcells)
  mats.filt.stemcells <- mats.filt.merge.mark[, cols.stemcells]
  print(dim(mats.filt.stemcells))
  # assertthat::assert_that(ncol(mats.filt.stemcells) > 0)
  
  if (ncol(mats.filt.unenriched) > 0 & ncol(mats.filt.linneg) > 0){
    mats.filt.unenrichedxlinneg <- cbind(mats.filt.unenriched, mats.filt.linneg)
  } else {
    mats.filt.unenrichedxlinneg <- Matrix(nrow = 0, ncol = 0)
  }
  
  if (ncol(mats.filt.unenriched) > 0 & ncol(mats.filt.stemcells) > 0){
    mats.filt.unenrichedxstemcells <- cbind(mats.filt.unenriched, mats.filt.stemcells)
  } else {
    mats.filt.unenrichedxstemcells <- Matrix(nrow = 0, ncol = 0)
  }
  
  # # unenriched and linneg
  # cols.unenrichedxlinneg <- grepl(paste0("^", jgrep.unenriched), colnames(mats.filt.merge.mark)) | grepl(paste0("^", jgrep.linneg), colnames(mats.filt.merge.mark))
  # mats.filt.unenrichedxlinneg <- mats.filt.merge.mark[, cols.unenrichedxlinneg]
  # print(dim(mats.filt.unenrichedxlinneg))
  # # assertthat::assert_that(ncol(mats.filt.unenrichedxlinneg) > 0)
  # 
  # 
  # # unenriched and stem cell
  # cols.unenrichedxstemcells <- grepl(paste0("^", jgrep.unenriched), colnames(mats.filt.merge.mark)) | grepl(paste0("^", jgrep.stemcells), colnames(mats.filt.merge.mark))
  # mats.filt.unenrichedxstemcells <- mats.filt.merge.mark[, cols.unenrichedxstemcells]
  # print(dim(mats.filt.unenrichedxstemcells))
  # # assertthat::assert_that(ncol(mats.filt.unenrichedxstemcells) > 0)
  
  
  print(paste("Saving mats for mark:", jmark))
  # write to outputs
  if (!file.exists(outf.all)){
    if (ncol(mats.filt.merge.mark) > 0){
      print("Dim for all")
      print(dim(mats.filt.merge.mark))
      saveRDS(mats.filt.merge.mark, file = outf.all)
    } else {
      print("Skipping, Empty mats.filt.merge.mark")
    }
  }
  if (!file.exists(outf.unenriched)){
    if (ncol(mats.filt.unenriched) > 0){
      print("Dim for unenriched")
      print(dim(mats.filt.unenriched))
      saveRDS(mats.filt.unenriched, file = outf.unenriched)
    } else {
      print("Skipping, Empty mats.filt.unenriched")
    }
  }
  if (!file.exists(outf.linneg)){
    if (ncol(mats.filt.linneg) > 0){
      print("Dim for linneg")
      print(dim(mats.filt.linneg))
      saveRDS(mats.filt.linneg, file = outf.linneg)
    } else {
      print("Skipping, Empty mats.filt.linneg")
    }
  }
  if (!file.exists(outf.stemcells)){
    if (ncol(mats.filt.stemcells) > 0){
      print("Dim for stemcells")
      print(dim(mats.filt.stemcells))
      saveRDS(mats.filt.stemcells, file = outf.stemcells)
    } else {
      print("Skipping, Empty mats.filt.stemcells")
    }
  }
  if (!file.exists(outf.unenrichedxlinneg)){
    if (ncol(mats.filt.stemcells) > 0){
      print("Dim for unenriched X linneg")
      print(dim(mats.filt.unenrichedxlinneg))
      saveRDS(mats.filt.unenrichedxlinneg, file = outf.unenrichedxlinneg)
    } else {
      print("Skipping, Empty outf.unenrichedxlinneg")
    }
  }
  if (!file.exists(outf.unenrichedxstemcells)){
    if (ncol(mats.filt.unenrichedxstemcells) > 0){
      print("Dim for unenriched X stemcells")
      print(dim(mats.filt.unenrichedxstemcells))
      saveRDS(mats.filt.unenrichedxstemcells, file = outf.unenrichedxstemcells)
    } else {
      print("Skipping, Empty mats.filt.unenrichedxstemcells")
    }
  }
}






