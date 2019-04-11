# Jake Yeung
# summarize_activities_all.R
# 2019-04-10
# DESCRIPTION
# 
#     Summarize activities 
# 
# FOR HELP
# 
#     Rscript summarize_activities_all.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-10
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

library(dplyr)
library(tidyr)
library(data.table)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('indir', metavar='INDIR',
                                            help='MARA output')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='output activities in long format')
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

inf.act <- file.path(args$indir, "Activities.gz")
inf.colname <- file.path(args$indir, "Colnames.gz")

print("Reading dat")
dat <- read.table(gzfile(inf.act), header=FALSE)
print("Reading cnames")
cnames <- unlist(read.table(gzfile(inf.colname), header=FALSE), use.names=FALSE)

colnames(dat) <- c("motif", cnames)

print("Gathering columns")
act.long <- tidyr::gather(dat, key = cell, value = activity, -motif)

print("Writing long dat")
fwrite(act.long, file = args$outfile)

system(paste("gzip", args$outfile))
