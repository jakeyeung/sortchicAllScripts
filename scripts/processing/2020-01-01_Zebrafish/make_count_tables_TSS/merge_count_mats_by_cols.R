# Jake Yeung
# merge_count_mats_by_cols.R
# 2019-11-08
# DESCRIPTION
# 
#     Merge counts by col
# 
# FOR HELP
# 
#     Rscript merge_count_mats_by_cols.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-11-08
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(Matrix)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('--infile', metavar='INFILE', nargs="+", required=TRUE, 
                                            help='List of input .RData to merge')
parser$add_argument('--outfile', metavar='OUTFILE', required=TRUE, 
                                            help='Output .RData')
# parser$add_argument('--exclude', metavar='EXCLUDE', required=FALSE,
#                                             help='Exclude infile starting with this string')
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

mats <- lapply(args$infile, function(inf){
  load(inf)  # count.dat
  return(count.dat$counts)
})
# get common rows
all.rows <- lapply(mats, function(mat){
  return(rownames(mat))
})

common.rows <- Reduce(intersect, all.rows)

mats <- lapply(mats, function(mat){
  mat[common.rows, ]
})

mats.merged <- do.call(cbind, mats)

count.dat <- list()
count.dat$counts <- mats.merged

save(count.dat, file = args$outfile)
