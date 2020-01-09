# Jake Yeung
# Date of Creation: 2019-12-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/0-get_good_cells_good_bins.noscrapping.R
# no scraping



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Load data ---------------------------------------------------------------

inf.meta <- "/home/jyeung/hpc/intestinal_scchic/metadata_all/OUD_to_filenames.txt"
meta <- fread(inf.meta, header = FALSE, col.names = c("Run", "Fastq")) %>%
  rowwise() %>%
  mutate(experi = strsplit(Fastq, "_")[[1]][[1]]) %>%
  dplyr::select(-Fastq)
meta <- meta[!duplicated(meta), ]

indir <- "/home/jyeung/hpc/intestinal_scchic/raw_data/HVG-intestines.tagged_bams.all_merged/countTables"

# jmark <- "k4me1"
jdate <- "2019-12-20"

# get LH counts
jname <- "NoScraping"
jgrp <- paste0("HVG_NoScraping_Bl6_.*.LHcounts.csv")
(infs.lh <- list.files(path = indir, pattern = jgrp, all.files = TRUE, full.names = TRUE))

# assertthat::assert_that(file.exists(infs.lh))
sapply(infs.lh, function(inf) assertthat::assert_that(file.exists(inf)))

jmarks <- sapply(infs.lh, function(x) strsplit(basename(x), split = "_")[[1]][[4]])
names(infs.lh) <- jmarks




# Constants ---------------------------------------------------------------


counts.cutoff.low <- 500  # H3K4me3
counts.cutoff.high <- 1000   # other marks 
TA.cutoff <- 0.5
counts.cutoff.lst <- vector(mode = "list", length = length(jmarks))
names(counts.cutoff.lst) <- jmarks



outmain <- "/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/quality_controls"
outpdf <- file.path(outmain, paste0("quality_controls_", jname, "_all_marks.", Sys.Date(), ".pdf"))

# Plot data ---------------------------------------------------------------

dats.rz <- ReadLH.SummarizeTA(infs.lh, remove.nones = FALSE, na.to.zero = TRUE, bind.rows = FALSE)

# label mark
dats.rz <- lapply(jmarks, function(jmark){
  jtmp <- dats.rz[[jmark]]
  jtmp$mark <- jmark
  return(jtmp)
}) %>%
  bind_rows() %>%
  rowwise() %>% 
  mutate(experi = ClipLast(samp, jsep = "_"))

# label runs 
dats.rz <- left_join(dats.rz, meta) %>%
  rowwise() %>%
  mutate(experiRun = ifelse(startsWith(experi, "HVG-"), paste(Run, experi, sep = "-"), experi))

# plot good cels 
m.all <- ggplot(dats.rz, aes(x = log10(total.count), y = TA.frac)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

m.dens <- ggplot(dats.rz, aes(x = log10(total.count), fill = mark)) + geom_density() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1)

# split by Lgr5? 
for (jmark in jmarks){
  print(jmark)
  x <- table(subset(dats.rz, mark == jmark)$experiRun)
  print(x)
}

dat.sum <- dats.rz %>%
  group_by(experiRun) %>%
  summarise(ncells = length(samp))

dat.sum2 <- dats.rz %>%
  group_by(experi) %>%
  summarise(ncells = length(samp)) %>%
  arrange(desc(ncells))

hist(dat.sum$ncells)

# add empty wells
empty.wells <- GetEmptyWells()

dats.rz$cellindx <- sapply(dats.rz$samp, function(x) paste("cell", strsplit(x, "_")[[1]][[2]], sep = ""))

dats.rz$is.empty <- dats.rz$cellindx %in% empty.wells

marks.with.low.counts <- c("k4me3", "k36me3")
for (jmark in jmarks){
  counts.cutoff.lst[[jmark]] <- ifelse(jmark %in% marks.with.low.counts, counts.cutoff.low, counts.cutoff.high)
}

dats.rz <- dats.rz %>%
  rowwise() %>%
  mutate(is.good = total.count > counts.cutoff.lst[[mark]] & TA.frac > TA.cutoff & !is.empty)


pdf(file = outpdf, useDingbats = FALSE)
print(m.all)
print(m.dens)
for (jmark in jmarks){
  m <- ggplot(dats.rz %>% filter(mark == jmark), aes(x = log10(total.count), y = TA.frac, color = is.good, size = is.empty, shape = is.empty)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~interaction(experiRun)) + ggtitle(jmark)
  print(m)
}
dev.off()

table(subset(dats.rz, mark == "k27me3")$experiRun)


# Define cells.keep -------------------------------------------------------

cells.keep <- subset(dats.rz, is.good)$samp

