# Jake Yeung
# Date of Creation: 2020-02-01
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/7-merge_new_H3K4me1_linneg_with_H3K4me1.R
# New experiment came in for H3K4me1 Linneg. Merge with other data and write file with
# removed plates and all plates. 
# Filter low variance cells beforehand. 
# Use bins from orignial 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)



# Settings ----------------------------------------------------------------

mincounts <- 10^3
minTAfrac <- 0.5
min.var <- 0.3

# Load older dat to be added ----------------------------------------------


# inf.mat.orig <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM.2020-01-31/BM_H3K4me1.varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.rds"
inf.mat.orig <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM.2020-01-31/BM_H3K4me1.varcutoff_0.3.platesRemoved.SmoothBinSize_1000.2020-02-04.Unenriched.rds"
count.mat.orig <- readRDS(inf.mat.orig)

bnameout <- strsplit(basename(inf.mat.orig), ".rds")[[1]][[1]]


# Load new count tables ---------------------------------------------------

# filter out low cells first

inf.lh <- "/home/jyeung/hpc/scChiC/raw_data/ZellerRawData_B6_linneg_K4me1.tagged_bams_mergedbymarks/LHcounts/H3K4me1-BM_Linneg_SC-merged.2020-01-31.tagged.LH_counts.csv"
dat.lh <- ReadLH.SummarizeTA(inf.lh) %>%
  rowwise() %>%
  mutate(plate = ClipLast(samp, "_"))

ggplot(dat.lh, aes(x = total.count, y = TA.frac)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() 

ggplot(dat.lh, aes(x = total.count, y = TA.frac)) + geom_point(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10()  + facet_wrap(~plate)

# filter out new plates only 

empty.wells <- GetEmptyWells()
dat.lh.filt <- subset(dat.lh, grepl("BM-lin-H3K4me1", samp)) %>%
  rowwise() %>%
  mutate(cellindx = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
         is.empty = cellindx %in% empty.wells,
         is.good = total.count >= mincounts & TA.frac >= minTAfrac)

ggplot(dat.lh.filt, aes(x = total.count, y = TA.frac, size = is.empty, color = is.good)) + geom_point(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_x_log10()  + facet_wrap(~plate)

cells.filt <- subset(dat.lh.filt, !is.empty & is.good)$samp
bins.keep <- rownames(count.mat.orig)

# Filt low var cells  -----------------------------------------------------

binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * mergesize



inf.mat <- "/home/jyeung/hpc/scChiC/raw_data/ZellerRawData_B6_linneg_K4me1.tagged_bams_mergedbymarks/countTables/H3K4me1-BM_Linneg_SC-merged.2020-01-31.tagged.countTable.csv"
mat.new <- ReadMatSlideWinFormat(inf.mat, as.sparse = TRUE)

mat.new.filt <- as.matrix(mat.new)[bins.keep, cells.filt]

mat.merge <- cbind(count.mat.orig, mat.new.filt)

dat.var.raw <- CalculateVarRaw(mat.new.filt, merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
ggplot(dat.var.raw, aes(x = ncuts, y = ncuts.var)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + geom_hline(yintercept = min.var)

dat.var.raw.merge <- CalculateVarRaw(mat.merge, merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

ggplot(dat.var.raw.merge, aes(x = ncuts, y = ncuts.var)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + facet_wrap(~plate) + scale_y_log10()


dat.lh.filt.merge <- left_join(dat.lh.filt, dat.var.raw, by = c("samp" = "cell"))

dat.lh.filt.merge <- dat.lh.filt.merge %>%
  rowwise() %>%
  mutate(is.lowvar = ncuts.var <= min.var) %>%
  mutate(IsGood.NotEmpty.NotLowVar = is.good & !is.empty & !is.lowvar)


cells.keep.final <- subset(dat.lh.filt.merge, IsGood.NotEmpty.NotLowVar)$samp

print(paste("Final ncells kept:", length(cells.keep.final)))



# Save to output ----------------------------------------------------------

count.mat.merged.filt <- cbind(count.mat.orig, mat.new.filt[, cells.keep.final])


outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM_includeNewLinnegH3K4me1.2020-01-31"
outf <- file.path(outdir, paste0(bnameout, ".WithLinneg.rds"))
saveRDS(count.mat.merged.filt, file = outf)


