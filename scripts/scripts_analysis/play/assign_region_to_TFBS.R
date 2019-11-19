# Jake Yeung
# Date of Creation: 2019-03-07
# File: ~/projects/scchic/scripts/scripts_analysis/play/assign_region_to_TFBS.R
# Assign region to TFBS so we can remove promoter regions. 


library(data.table)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")
source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Read section of table ---------------------------------------------------

inf <- "/private/tmp/test.bed"

dat <- data.table::fread(inf)

colnames(dat) <- c("chromo", "start", "end", "motif", "count", "gene", "dist", "peak")


# Make genomic regions and assign reads -----------------------------------

jsub <- subset(dat, motif == "Ahr")
regions <- data.frame(seqnames = sapply(jsub$peak, GetChromo),
                      start = jsub$start,
                      end = jsub$end, 
                      peak = jsub$peak,
                      motif = jsub$motif,
                      count = jsub$count,
                      gene = jsub$gene,
                      stringsAsFactors = FALSE)
# rownames(regions) <- make.names(jsub$peak, unique = TRUE, names = FALSE)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
# rownames(regions.range) <- make.names(regions$peak, unique = TRUE, names = FALSE)
regions.range$motif <- regions$motif
regions.range$peak <- regions$peak
regions.range$count <- regions$count
regions.range$gene <- regions$gene
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db', assignGenomicAnnotation = FALSE, verbose = FALSE))
data.table::fwrite(subset(regions.annotated, select = c(seqnames, start, end, motif, count, gene, distanceToTSS, peak)), 
                   file = "/tmp/annot.bed", quote = FALSE, sep = "\t", col.names = FALSE)
# regions.annotated$region_coord <- names(regions.range)

# top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))
