# Jake Yeung
# Date of Creation: 2019-06-17
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/1-quality_control_get_good_cells.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

source("scripts/Rfunctions/QCFunctions.R")



# Load data ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

prefix <- "_PZ-Bl6-BM-Linneg"
# indir <- "/Users/yeung/data/scchic/from_cluster/TA_count_frequencies_bugfixed_uniqname"
indir <- paste0("/Users/yeung/data/scchic/from_cluster/count_mat_from_bam", prefix)
infs <- list.files(indir, pattern = "*.RZcounts.csv", full.names = TRUE)

system.time(
  dat.long <- lapply(infs, ReadDinucLinneg) %>%
    bind_rows()
)



# Count cells -------------------------------------------------------------

dat.counts.sum <- dat.long %>%
  group_by(samp, mark) %>%
  summarise(count = sum(count))

ggplot(dat.counts.sum, aes(x = log10(count))) + geom_histogram() + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Now we count TA frequency  ----------------------------------------------

dat.dinuc.sum <- dat.long %>%
  group_by(samp, mark) %>%
  mutate(dinuc.frac = count / sum(count)) %>%
  filter(dinuc == "TA")


# Now plot TA versus counts -----------------------------------------------

empty.wells <- GetEmptyWells()
dat.sum.merged <- left_join(dat.counts.sum, dat.dinuc.sum %>% dplyr::select(samp, mark, dinuc, dinuc.frac, cell))

dat.sum.merged$is.empty <- sapply(dat.sum.merged$cell, function(x) ifelse(x %in% empty.wells, TRUE, FALSE))

ggplot(dat.sum.merged, aes(x = log10(count), y = dinuc.frac, color = is.empty)) + geom_point() + 
  facet_wrap(~mark, nrow = 1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# get counts per cell
dat.counts.sum <- dat.long %>%
  group_by(samp, mark) %>%
  summarise(count = sum(count))
dat.counts.sum$mark <- factor(dat.counts.sum$mark, levels = jmarks)

ggplot(dat.counts.sum, aes(x = log10(count))) + geom_histogram(bins = 60) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# plot fraction of TA
dat.frac.sum <- dat.long %>%
  group_by(samp, mark) %>%
  mutate(dinuc.frac = count / sum(count)) %>%
  filter(dinuc == "TA") %>%
  dplyr::select(dinuc, samp, mark, dinuc.frac, cell)

# get empty wells 
empty.wells <- GetEmptyWells(indx = 0)

dat.sum.merged <- left_join(dat.counts.sum, dat.frac.sum)

dat.sum.merged$is.empty <- sapply(dat.sum.merged$cell, function(x) ifelse(x %in% empty.wells, TRUE, FALSE))

ggplot(dat.sum.merged, aes(y = dinuc.frac, x = log10(count), color = is.empty)) + geom_point() + 
  facet_wrap(~mark, nrow = 1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Find good cells and bad cells  ------------------------------------------

# p <- 0.95  # cell counts
# p.stringent <- 0.9  # cell counts
# pfrac <- 0.95  # TA fraction
# be more stringent with H3K4me3?
TAcutoff <- 0.7
countcutoff <- 3

dat.sum.merged.goodbad <- dat.sum.merged %>%
  mutate(cellsum.log = log10(count)) %>%
  group_by(mark) %>%
  mutate(good.cellsize = ifelse(cellsum.log > countcutoff, TRUE, FALSE),
         good.frac = ifelse(dinuc.frac > TAcutoff, TRUE, FALSE),
         good.cell = as.logical(good.cellsize * good.frac))

# dat.sum.merged.relaxed <- dat.sum.merged %>%
#   mutate(cellsum.log = log10(count)) %>%
#   group_by(mark) %>%
#   mutate(good.cellsize = ifelse(cellsum.log > quantile(cellsum.log, probs = (1-p)), TRUE, FALSE),
#          good.frac = ifelse(dinuc.frac > quantile(dinuc.frac, probs = (1-pfrac)), TRUE, FALSE),
#          good.cell = as.logical(good.cellsize * good.frac)) %>%
#   dplyr::filter(mark %in% c("H3K27me3", "H3K9me3"))
# 
# dat.sum.merged.stringent <- dat.sum.merged %>%
#   mutate(cellsum.log = log10(count)) %>%
#   group_by(mark) %>%
#   mutate(good.cellsize = ifelse(cellsum.log > quantile(cellsum.log, probs = (1-p.stringent)), TRUE, FALSE),
#          good.frac = ifelse(dinuc.frac > quantile(dinuc.frac, probs = (1-pfrac)), TRUE, FALSE),
#          good.cell = as.logical(good.cellsize * good.frac)) %>%
#   dplyr::filter(mark %in% c("H3K4me1", "H3K4me3"))
# dat.sum.merged.goodbad <- bind_rows(dat.sum.merged.relaxed, dat.sum.merged.stringent)

dat.sum.merged.goodbad$mark <- factor(as.character(dat.sum.merged.goodbad$mark), levels = jmarks)

ggplot(dat.sum.merged.goodbad, aes(x = cellsum.log, y = dinuc.frac, color = good.cell, shape = is.empty)) + geom_point(size = 2) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Print good cells --------------------------------------------------------

good.cells <- data.frame(cell = subset(dat.sum.merged.goodbad, good.cell)$samp)
# outtxt <- paste0("~/data/scchic/quality_control", prefix, "/good_cells.pcounts.", (1-p), ".pfrac.", (1-pfrac), ".pcountsStringent.", (1-p.stringent), ".txt")
outtxt <- paste0("~/data/scchic/quality_control", prefix, "/good_cells.TAcutoff.", TAcutoff, ".countcutoff.", countcutoff, ".txt")
if (!file.exists(outtxt)){
  fwrite(good.cells, file = outtxt, col.names = FALSE)
}





