# Jake Yeung
# Date of Creation: 2018-12-15
# File: ~/projects/scChiC/scripts/quality_controls/summarize_barcodes.R
# Summarize barcodes

args <- commandArgs(trailingOnly=TRUE)

print(args)
barcodes.dir <- args[[1]]
count.thres <- as.numeric(args[[2]])
outdir <- args[[3]]

if (is.na(count.thres)){
  stop(paste("Count thres must be numeric", count.thres))
}

library(dplyr)
library(ggplot2)
library(stringr)

source("scripts/Rfunctions/GetMetaData.R")

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

print(paste("Current directory: ", getwd()))

# infiles <- list.files("data/barcode_summaries", pattern = "*bc_counts.txt", full.names = TRUE)
infiles <- list.files(barcodes.dir, pattern = "*bc_counts.txt", full.names = TRUE)

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

m1 <- ggplot(dat, aes(x = jrank, y = counts, color = biorep)) + 
  geom_point() + facet_grid(experi ~ chip) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

m2 <- ggplot(dat, aes(x = jrank, y = countsfrac, color = biorep)) + 
  geom_point() + facet_grid(experi ~ chip) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10()

pdf(file.path(outdir, "diagnostic_plots.pdf"), useDingbats=FALSE)
  print(m1)
  print(m2)
dev.off()

jchips=c("H3K27me3", "H3K4me1", "H3K4me3", "H3K9me3")
count.thres <- 0

# find cutoff
# jchip <- "H3K27me3"

for (jchip in jchips){
  # apply stringent threshold 
  jsub <- subset(dat, chip == jchip & counts > count.thres)
  dim(jsub)  # ~1033 samples at 10000 threshold
  
  
  # outdir <- "outputs_R/barcode_summaries"
  # dir.create(outdir)
  # write down summarizes files 
  write.table(dat, file = file.path(outdir, paste0("barcode_summaries_all.", jchip, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(jsub, file = file.path(outdir, paste0("barcode_summaries_filtered.", count.thres, ".chip.", jchip, ".txt")), quote = FALSE, sep = "\t", row.names = FALSE)
  
  # make a summary per fbase (easier downstream processing)
  jsub.by.bam <- split(jsub, jsub$fbase)
  
  lapply(jsub.by.bam, function(jsubsplit){
    write.table(jsubsplit, file = 
                  file.path(outdir, 
                            paste0("barcode_summary.", unique(jsubsplit$fbase), ".thres.", count.thres, ".chip.", jchip, ".txt")), quote = FALSE, sep = "\t", row.names=FALSE)
  })
  
}

