# Jake Yeung
# annotate_distance_on_bed.R
# 2019-03-07
# DESCRIPTION
# 
#     Annotate distance on bed 
# 
# FOR HELP
# 
#     Rscript annotate_distance_on_bed.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-03-07
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time()

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='inbed')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='outbed')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

library(data.table)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
source("scripts/Rfunctions/Aux.R")  # set WD properly

# Do stuff
dat <- data.table::fread(args$infile)
print("Finished reading file...")
colnames(dat) <- c("chromo", "start", "end", "motif", "count", "gene", "dist", "peak")

print("Preparing genomic regions...")
regions <- data.frame(seqnames = sapply(dat$peak, GetChromo),
                      start = dat$start,
                      end = dat$end, 
                      peak = dat$peak,
                      motif = dat$motif,
                      count = dat$count,
                      gene = dat$gene,
                      stringsAsFactors = FALSE)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.range$motif <- regions$motif
regions.range$peak <- regions$peak
regions.range$count <- regions$count
regions.range$gene <- regions$gene
print("Annotating regions...")
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db', assignGenomicAnnotation = FALSE, verbose = FALSE))
data.table::fwrite(subset(regions.annotated, select = c(seqnames, start, end, motif, count, gene, distanceToTSS, peak)), 
                   file = args$outfile, quote = FALSE, sep = "\t", col.names = FALSE)

print(Sys.time() - jstart)
