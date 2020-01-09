# Jake Yeung
# filter_mat_by_good_cells.R
# 2019-11-08
# DESCRIPTION
# 
#     Filter sparse mat by good cells
# 
# FOR HELP
# 
#     Rscript filter_mat_by_good_cells.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-11-08
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

library(Matrix)
library(data.table)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='.rds')
parser$add_argument('goodcells', metavar='Column of good cells',
                                            help='.csv')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='.RData named count.dat$counts filtered')
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

# get good cells 
dat <- data.table::fread(args$goodcells, header=FALSE)  # V1

print(dat)
cells <- dat$V1

assertthat::assert_that(length(cells) > 0)

print(cells)

count.mat <- readRDS(args$infile)

cols.keep <- which(colnames(count.mat) %in% cells)

print(colnames(count.mat)[cols.keep])

count.dat <- list()
count.dat$counts <- count.mat[, cols.keep]
if (ncol(count.dat$counts) > 0){
  save(count.dat, file = args$outfile)
} else {
  print("no good cells matched, doing nothing")
  stop("No good cells matched, exiting")
}








