# Jake Yeung
# summarize_bedfiles.R
# 2021-06-15
# DESCRIPTION
# 
#     Summarize bed
# 
# FOR HELP
# 
#     Rscript summarize_bedfiles.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2021-06-15
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input bed.gz')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output bed summary')
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

cnames <- c("chromo", "start", "end", "dupcounts", "bc")
dat.tmp <- fread(args$infile, header = FALSE, col.names = cnames)

fname <- basename(args$infile)
fbase <- strsplit(fname, split = "\\.")[[1]][[1]]

dat.tmp$cellname <- paste(fbase, dat.tmp$bc, sep = "_") 
dat.tmp$uniquecounts <- 1

dat.tmp.sum <- dat.tmp %>%
  group_by(cellname) %>%
  summarise(dupcounts = sum(dupcounts),
            uniquecounts = sum(uniquecounts)) %>%
  arrange(desc(uniquecounts))
dat.tmp.sum$fbase <- fbase

# write output
fwrite(dat.tmp.sum, file = args$outfile, sep = "\t", col.names = TRUE)
