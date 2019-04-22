# Jake Yeung
# make_bed_from_count_mat.R
# 2019-04-19
# DESCRIPTION
# 
#     Make bed from count mat
# 
# FOR HELP
# 
#     Rscript make_bed_from_count_mat.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-19
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(Matrix)

source("scripts/Rfunctions/Aux.R")

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='RData',
                                            help='.RData')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='bed file')
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

load(args$infile, v=T)

rnames <- sapply(rownames(count.dat$counts), function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE)
print(head(rnames))
chromos <- sapply(rnames, GetChromo)
starts <- sapply(rnames, GetStart)
ends <- sapply(rnames, GetEnd)

print(head(chromos))
dat <- data.frame(chromo = chromos,
                  start = starts,
                  end = ends)

print(head(dat))

data.table::fwrite(dat, file = args$outfile, sep = "\t", col.names = FALSE, row.names = FALSE)

