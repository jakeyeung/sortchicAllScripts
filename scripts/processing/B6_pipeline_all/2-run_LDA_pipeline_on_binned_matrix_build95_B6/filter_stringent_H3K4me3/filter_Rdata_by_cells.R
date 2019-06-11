# Jake Yeung
# filter_Rdata_by_cells.R
# 2019-06-03
# DESCRIPTION
# 
#     Filter count matrix by cells 
# 
# FOR HELP
# 
#     Rscript filter_Rdata_by_cells.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-03
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

library(Matrix)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='RData')
parser$add_argument('cellfile', metavar='CELLFILE',
                                            help='List of cells. Includes header')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='RData filtered by cells')
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

load(args$infile, v=T)  # count.dat$counts
cells <- read.table(args$cellfile, header=TRUE)

cols.keep <- which(colnames(count.dat$counts) %in% cells$cell)
count.dat$counts <- count.dat$counts[, cols.keep]
save(count.dat, file = args$outfile)
