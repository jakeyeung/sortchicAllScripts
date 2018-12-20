# Jake Yeung
# Date of Creation: 2018-12-15
# File: ~/projects/scChiC/scripts/quality_controls/summarize_barcodes.R
# Summarize barcodes

library(dplyr)
library(ggplot2)
library(stringr)

GetTissue <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> BM
  x <- strsplit(inf, split = "-")[[1]][[2]]
  return(x)
}
GetChip <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H3K7me3
  x <- strsplit(inf, split = "-")[[1]][[4]]
  return(x)
}
GetBioRep <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> m1 -> 1
  x <- strsplit(inf, split = "-")[[1]][[3]]
  x <- strsplit(x, split = "m")[[1]][[2]]
  return(x)
}
GetTechRep <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> S14-> 14
  x <- stringr::str_match(inf, "_S[1-9]*.")[1, 1]
  # remove first and last characters
  x <- substr(x, start = 2, stop = nchar(x) - 1)
  return(x)
}
GetExperiment <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H2GV2BGX9
  x <- strsplit(inf, split = "_")[[1]][[2]]
  return(x)
}

AddMetaToDat <- function(dattmp, inf){
  # add meta data to columns of dat
  dattmp$tissue <- GetTissue(inf)
  dattmp$chip <- GetChip(inf)
  dattmp$biorep <- GetBioRep(inf)
  dattmp$techrep <- GetTechRep(inf)
  dattmp$experi <- GetExperiment(inf)
  dattmp$fbase <- strsplit(inf, "\\.")[[1]][[1]]
  return(dattmp)
}

# Summarize barcodes ------------------------------------------------------

infiles <- list.files("data/barcode_summaries", pattern = "*bc_counts.txt", full.names = TRUE)

dat <- lapply(infiles, function(inf){
  dattmp <- read.table(inf, header = FALSE, col.names = c("counts", "barcode"))
  dattmp <- AddMetaToDat(dattmp, basename(inf))
  return(dattmp)
}) %>%
  dplyr::bind_rows(.) %>%
  group_by(fbase) %>%
  arrange(desc(counts)) %>%
  mutate(jrank = seq(length(counts)),
         countsfrac = (counts) / sum(counts))

ggplot(dat, aes(x = jrank, y = counts, color = biorep)) + 
  geom_point() + facet_grid(experi ~ chip) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

ggplot(dat, aes(x = jrank, y = countsfrac, color = biorep)) + 
  geom_point() + facet_grid(experi ~ chip) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

# Let’s pick H3K27me3 -----------------------------------------------------

# find cutoff
jchip <- "H3K27me3"
count.thres <- 10000
# apply stringent threshold 
jsub <- subset(dat, chip == jchip & counts > count.thres)
dim(jsub)  # ~1033 samples at 10000 threshold


outdir <- "outputs_R/barcode_summaries"
dir.create(outdir)
# write down summarizes files 
write.table(dat, file = file.path(outdir, "barcode_summaries_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(jsub, file = file.path(outdir, paste0("barcode_summaries_filtered.", count.thres, ".chip.", jchip, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)

# make a summary per fbase (easier downstream processing)
jsub.by.bam <- split(jsub, jsub$fbase)

lapply(jsub.by.bam, function(jsubsplit){
  write.table(jsubsplit, file = 
                file.path(outdir, 
                          paste0("barcode_summary.", unique(jsubsplit$fbase), ".thres.", count.thres, ".chip.", jchip, ".txt")), quote = FALSE, sep = "\t", row.names=FALSE)
})

