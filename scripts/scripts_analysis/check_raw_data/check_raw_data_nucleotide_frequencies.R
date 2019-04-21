# Jake Yeung
# Date of Creation: 2019-04-16
# File: ~/projects/scchic/scripts/scripts_analysis/check_raw_data/check_raw_data_nucleotide_frequencies.R
# Check raw nucleotide frequencies

rm(list=ls())

library(dplyr)
library(ggplot2)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")

# Functions ---------------------------------------------------------------

ParseFilename <- function(fname){
  # make into H3K4me1.BM.m1.S9.AH3VGVBGX9.CAGAATAT so we can use MARA loading scripts to convert names
  # fname format: PZ-BM-m2-H3K27me3-2_H2GV2BGX9_S18_L002_R1_001.summary.out
  fbase <- strsplit(fname, split = "\\.")[[1]][[1]]
  mark <- strsplit(fbase, split = "-")[[1]][[4]]
  mouse <- strsplit(fbase, split = "-")[[1]][[3]]
  tiss <- strsplit(fbase, split = "-")[[1]][[2]]
  techrep <- strsplit(fbase, split = "_")[[1]][[3]]
  techrepname <- strsplit(fbase, split = "_")[[1]][[2]]
  cellname <- paste(mark, tiss, mouse, techrep, techrepname, sep = ".")  # no cell barcode need to add it yourself 
  return(cellname)
}

ReadAndParseFile <- function(inf){
  dat <- data.table::fread(inf, header = FALSE, stringsAsFactors = FALSE, col.names = c("cellname", "bc", "dinuc"))
  return(dat)
}

FilterAndSwitchColnames <- function(dat.sum){
  bcs.uniq <- unique(dat.sum$bc)
  bcs.keep <- sapply(bcs.uniq, function(x){
    ifelse(!is.null(cellhash[[x]]), cellhash[[x]], NA)
  })
  bcs.keep <- names(bcs.keep[!is.na(bcs.keep)])
  dat.sum <- subset(dat.sum, bc %in% bcs.keep)
  dat.sum$cname.new <- SwitchColnames(dat.sum$cellname)
  return(dat.sum)
}



# Load files --------------------------------------------------------------


indir <- "/Users/yeung/data/scchic/from_cluster/count_summaries/parsed2"
indir2 <- "/Users/yeung/data/scchic/from_cluster/count_summaries_VAN2979"

infs <- list.files(indir, pattern = "*.parsed.out", full.names = TRUE)
infs2 <- list.files(indir2, pattern = "*.parsed.out", full.names = TRUE)

dat <- lapply(infs, function(inf){
  return(ReadAndParseFile(inf))
}) %>% 
  bind_rows()

dat2 <- lapply(infs2, function(inf){
  return(ReadAndParseFile(inf))
}) %>% 
  bind_rows()
  
dat.sum <- dat %>%
  group_by(cellname, bc, dinuc) %>%
  summarise(dinuc.freq = length(dinuc))
dat.sum <- dat.sum %>%
  group_by(cellname) %>%
  mutate(dinuc.frac = dinuc.freq / sum(dinuc.freq))

dat.sum2 <- dat2 %>%
  group_by(cellname, bc, dinuc) %>%
  summarise(dinuc.freq = length(dinuc))
dat.sum2 <- dat.sum2 %>%
  group_by(cellname) %>%
  mutate(dinuc.frac = dinuc.freq / sum(dinuc.freq))



# Add proper cell labels --------------------------------------------------

cellhash <- GetCellHashBC()



dat.sum2 <- FilterAndSwitchColnames(dat.sum2)

# bcs.uniq <- unique(dat.sum$bc)
# bcs.keep <- sapply(bcs.uniq, function(x){
#   ifelse(!is.null(cellhash[[x]]), cellhash[[x]], NA)
# })
# bcs.keep <- names(bcs.keep[!is.na(bcs.keep)])
# dat.sum <- subset(dat.sum, bc %in% bcs.keep)
# dat.sum$cname.new <- SwitchColnames(dat.sum$cellname)

dat.sum$mark <- sapply(dat.sum$cellname, function(x) strsplit(x, "\\.")[[1]][[1]])
dat.sum2$mark <- sapply(dat.sum2$cellname, function(x) strsplit(x, "\\.")[[1]][[1]])

# Summarize dinucloetide frequencies across cells and marks ---------------

# dat.filt <- subset(dat.sum, grepl("^TA$", dinuc))
dat.filt <- subset(dat.sum, dinuc == "TA")
dat.filt2 <- subset(dat.sum2, dinuc == "TA")

dat.filt.merged <- bind_rows(dat.filt, dat.filt2)

ggplot(dat.filt, aes(x = dinuc.frac)) + geom_histogram(aes(y=..count../(sum(..count..)))) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Fraction of TA")

ggplot(dat.filt.merged, aes(x = dinuc.frac)) + geom_histogram(aes(y=..count../(sum(..count..)))) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Fraction of TA") + facet_wrap(~mark)

ggplot(dat.filt.merged, aes(x = dinuc.frac)) + geom_histogram() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Fraction of TA") + facet_wrap(~mark) + geom_vline(data = dat.filt.merged %>% group_by(mark) %>% summarise(dinuc.frac.mean = median(dinuc.frac)), aes(xintercept = dinuc.frac.mean)) + 
  ggtitle("Vertical line is median")

# Load UMAP coordinates and project ---------------------------------------

load("/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-15.RData", v=T)

dat.trajs.long <- left_join(dat.trajs.long, dat.filt.merged %>% ungroup() %>% dplyr::select(cname.new, dinuc.frac), by = c("cell" = "cname.new"))

ggplot(dat.trajs.long, aes(x = dinuc.frac)) + geom_histogram() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Fraction of Total Cuts Starting with 'TA'") + facet_wrap(~mark)

# ggplot(dat.trajs.long, aes(x = dinuc.frac)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   xlab("Fraction of Total Cuts Starting with 'TA'") + facet_wrap(~mark)

ggplot(dat.trajs.long, aes(x = X1, y = X2, color = dinuc.frac)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + scale_color_gradient2(low = "gray85", mid = "gray50", high = muted("darkblue"), midpoint = 0.5)
