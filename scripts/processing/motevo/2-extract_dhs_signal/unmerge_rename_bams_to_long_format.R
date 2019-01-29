# Jake Yeung
# rename_bams_to_long_format.R
# Merge bams and rename into sensible long format names (ZT0 to 24 is NSD, 24 and on is SD)
# Same as nmerge_rename_bams_to_long_format.R" but no merging, just renaming
# 2016-07-16

library(stringr)
library(dplyr)
library(hash)

args <- commandArgs(trailingOnly=TRUE)
bamdir <- args[[1]]
outdir <- args[[2]]

print(paste("Bamdir:", bamdir))
print(paste("Outdir:", outdir))

if(!dir.exists(bamdir)) stop(paste("Bamdir does not exist:", bamdir))
dir.create(outdir)

# source
source("/home/jyeung/projects/sleep_deprivation/handle_bams/utils/HandleRawDataFunctions.R")
source("/home/jyeung/projects/sleep_deprivation/handle_bams/utils/StringFunctions.R")

# get metadata
dat.meta <- ReadProcessMetadata("/archive/epfl/upnae/jyeung/sleep_deprivation/SDrecovery_design_formatted.txt", keep.original.samps = TRUE)
# format metadata to long
dat.sd.shift <- subset(dat.meta, SD_NSD == "SD") %>%
  group_by(ZT) %>%
  mutate(time.shift = time + 24)
dat.sd.shift$ZT <- dat.sd.shift$time.shift; dat.sd.shift$time.shift <- NULL

dat.meta.shift <- rbind(subset(dat.meta, SD_NSD == "NSD"), as.data.frame(dat.sd.shift))


# make ZT0 ZT24
dat.meta.shift$ZT <- sapply(dat.meta.shift$ZT, function(zt) ifelse(zt == 0, 24, zt))
# add ZT and add leading zeros to time
dat.meta.shift$ZT <- sapply(dat.meta.shift$ZT, function(zt) paste("ZT", sprintf("%02d", zt), sep = ""))

dat.meta.shift.sub <- subset(dat.meta.shift, select = c(ZT, SD_NSD, samp))
new.names <- apply(dat.meta.shift.sub, 1, function(row) paste(row, collapse = "_")[[1]])
dat.meta.shift$new.names <- new.names

print(dat.meta.shift)

bamlist <- list.files(bamdir, pattern = "*bam*")
for (b in bamlist){
  b.base <- strsplit(b, "\\.")[[1]][[1]]
  ext <- paste(strsplit(b, "\\.")[[1]][-1], collapse = ".")
  b.new <- subset(dat.meta.shift, SampleID.orig == b.base)$new.names
  old.path <- file.path(bamdir, b)
  new.path <- file.path(outdir, paste0(b.new, ".", ext))
  cmd <- paste("ln -s", old.path, new.path)
  system(cmd)
}
