# Jake Yeung
# filter_bins_by_chromosomes.R
# 2020-04-08
# DESCRIPTION
# 
#     Keep only a list of chromosomes
# 
# FOR HELP
# 
#     Rscript filter_bins_by_chromosomes.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-08
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)


suppressPackageStartupMessages(library("argparse"))
library(Matrix)
library(scchicFuncs)
library(dplyr)
library(data.table)
library(ggplot2)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='in .rds count matrix, bad cells filtered out')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='out .rds, keeping bins corresponding to good chromosomes')
parser$add_argument('-chromoskeep', metavar='chromosomes to keep', nargs = "+",  
                                            help='Chromosomes to keep, space delimited text')
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

if (endsWith(args$infile, ".rds")){
  mat <- readRDS(args$infile)
} else if (endsWith(args$infile, ".csv")){
  mat <- ReadMatSlideWinFormat(args$infile, as.sparse=TRUE, add.coord=FALSE, sort.rnames=TRUE)
} else {
  stop(args$infile, "must end with .rds or .csv")
}


# show all chromos
jchromos.vec <- sapply(rownames(mat), function(x) strsplit(x, ":")[[1]][[1]])

print("All chromos:")
print(sort(unique(jchromos.vec)))

# filter out only good chromos?
jchromos.filt.i <- which(jchromos.vec %in% args$chromoskeep)

assertthat::assert_that(length(jchromos.filt.i) > 0)

mat.filt <- mat[jchromos.filt.i, ]

print(paste("After filtering for chromosomes, we have", length(jchromos.filt.i), "bins left"))

print(paste("Dimensions of new matrix:"))
print(dim(mat.filt))

# check chromos after
jchromos.check <- sort(unique(sapply(rownames(mat.filt), function(x) strsplit(x, ":")[[1]][[1]])))

print("Checking chromosomes after filtering matrix:")
print(jchromos.check)

# write to output
saveRDS(mat.filt, args$outfile)











