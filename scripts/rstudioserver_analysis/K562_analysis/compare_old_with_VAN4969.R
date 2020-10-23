# Jake Yeung
# Date of Creation: 2020-08-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/K562_analysis/compare_old_with_VAN4969.R
# Project stuff


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load old tables ---------------------------------------------------------

indir.old <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data/ZellerRawDataK562.tagged_bams_mergedbymarks/countTablesAndRZr1only.NewCountFilters"))
infs.rz <- list.files(path = indir.old, pattern = ".*.RZ.*.csv", full.names = TRUE)
infs.counts <- list.files(path = indir.old, pattern = ".*.countTable.winsize_50000*.csv", full.names = TRUE)
names(infs.rz) <- sapply(infs.rz, function(x) strsplit(basename(x), "-")[[1]][[4]])
names(infs.counts) <- sapply(infs.counts, function(x) strsplit(basename(x), "-")[[1]][[4]])

rz.old <- lapply(infs.rz, function(inf.rz){
  # inf.rz <- file.path(indir.old, paste0("PZ-K562-G1-", jmark, "-Rep1Rep2merged.tagged.RZsummary.csv"))
  dat <- ReadLH.SummarizeTA(inf.rz)
})

rz.old <- lapply(jmarks, function(jmark){
  rz.old[[jmark]]$mark <- jmark
  return(rz.old[[jmark]])
}) %>%
  bind_rows() 

empty.wells <- GetEmptyWells(indx = 0)

rz.old <- rz.old %>%
  rowwise() %>%
  mutate(wellindx = strsplit(samp, "_")[[1]][[2]], 
         cell = paste("cell", wellindx, sep = ""),
         is.empty = cell %in% empty.wells)


# filter cells 

ggplot(rz.old, aes(x = log10(total.count), y = TA.frac, color = is.empty, size = is.empty)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~mark) + 
  geom_vline(xintercept = log10(500))

# set cutoffs
cutoffs.dat <- data.frame(mark = jmarks, cutoff = c(1000, 500, 1000, 1000), stringsAsFactors = FALSE)


rz.old.merge <- left_join(rz.old, cutoffs.dat) %>%
  rowwise() %>%
  mutate(is.good = total.count >= cutoff)

cells.keep <- subset(rz.old.merge, !is.empty & is.good)$samp


# Load new data tables ----------------------------------------------------


mats.new <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, ".rds"))
  readRDS(inf)
})


lapply(mats.new, dim)



# load mats
mats.old <- lapply(infs.counts, function(inf){
  dat <- ReadMatSlideWinFormat(inf, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = TRUE)
  cols.keep <- colnames(dat) %in% cells.keep
  # rows.keep <- rownames(dat) %in% coords.keep
  # dat[rows.keep, cols.keep]
  return(dat)
})
lapply(mats.old, dim)


coords.keep.new <- Reduce(f = intersect, x = lapply(mats.new, function(jmat){
  rownames(jmat)
}))

coords.keep.old <- Reduce(f = intersect, x = lapply(mats.old, function(jmat){
  rownames(jmat)
}))

coords.keep <- intersect(coords.keep.new, coords.keep.old)


# Filter mats -------------------------------------------------------------

mats.old.filt <- lapply(mats.old, function(jmat){
  rows.keep <- rownames(jmat) %in% coords.keep
  return(jmat[rows.keep, ])
})

mats.new.filt <- lapply(mats.new, function(jmat){
  rows.keep <- rownames(jmat) %in% coords.keep
  return(jmat[rows.keep, ])
})


# Write old tables output ready for LDA -----------------------------------

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cleaned_count_tables_for_lda_and_projections"
# 
# for (jmark in jmarks){
#   outf.old <- file.path(outdir, paste0("count_mats_old_binsize_50000_genomewide.", jmark, ".old.rds"))
#   outf.new <- file.path(outdir, paste0("count_mats_old_binsize_50000_genomewide.", jmark, ".new.rds"))
#   if (!file.exists(outf.old)){
#     saveRDS(mats.old.filt[[jmark]], file = outf.old)
#   } else {
#     print(paste(outf.old, "exists, skipping"))
#   }
#   if (!file.exists(outf.new)){
#     # saveRDS(mats.old.filt[[jmark]], file = outf.old)
#     saveRDS(mats.new.filt[[jmark]], file = outf.new)
#   } else {
#     print(paste(outf.old, "exists, skipping"))
#   }
# }



# Filter top 5000 ---------------------------------------------------------


# load top 5000 

inf.test <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.H3K4me1.G1_G2_S.topn_5000.rds"
assertthat::assert_that(file.exists(inf.test))

mat.test <- readRDS(inf.test)

rows.top5000 <- intersect(rownames(mat.test), coords.keep)

outdir.top5000 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cleaned_count_tables_for_lda_and_projections.top5000.intersect"
dir.create(outdir.top5000)

for (jmark in jmarks){
  outf.old <- file.path(outdir.top5000, paste0("count_mats_old_binsize_50000_genomewide.top5000.", jmark, ".old.rds"))
  outf.new <- file.path(outdir.top5000, paste0("count_mats_old_binsize_50000_genomewide.top5000.", jmark, ".new.rds"))
  
  mat.old.tmp <- mats.old.filt[[jmark]]
  mat.new.tmp <- mats.new.filt[[jmark]]
  rows.keep.old <- rownames(mat.old.tmp) %in% rows.top5000
  rows.keep.new <- rownames(mat.new.tmp) %in% rows.top5000
  
  
  mat.old.tmp.filt <- mat.old.tmp[rows.keep.old, ]
  mat.new.tmp.filt <- mat.new.tmp[rows.keep.new, ]
  
  # remove empty cells
  cells.keep.old <- which(colSums(mat.old.tmp.filt) > 10)
  cells.keep.new <- which(colSums(mat.new.tmp.filt) > 10)
  
  mat.old.tmp.filt <- mat.old.tmp.filt[, cells.keep.old]
  mat.new.tmp.filt <- mat.new.tmp.filt[, cells.keep.new]
  
  print(jmark)
  
  # print(dim(mat.old.tmp.filt))
  # print(dim(mat.new.tmp.filt))
  
  print(range(rowSums(mat.old.tmp.filt)))
  print(range(rowSums(mat.new.tmp.filt)))
  
  print(range(colSums(mat.old.tmp.filt)))
  print(range(colSums(mat.new.tmp.filt)))
  
  saveRDS(mat.old.tmp.filt, file = outf.old)
  saveRDS(mat.new.tmp.filt, file = outf.new)
}


