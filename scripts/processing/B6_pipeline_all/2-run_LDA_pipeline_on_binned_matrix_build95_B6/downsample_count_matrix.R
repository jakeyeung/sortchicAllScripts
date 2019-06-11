# Jake Yeung
# 4_downsample_count_matrix.R
# 2019-06-03
# DESCRIPTION
# 
#     Downsample a matrix
# 
# FOR HELP
# 
#     Rscript downsample_count_matrix.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-03
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jtime <- Sys.time() 

library(DropletUtils)  # https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#from-the-count-matrix
library(Matrix)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Path to RData')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Path to output')
parser$add_argument('counts_final', metavar='counts_final', type='integer')
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

cellsizes <- colSums(count.dat$counts)
prop.vec <- args$counts_final / cellsizes
print(prop.vec)
cells.keep <- which(prop.vec <= 1)

print("Dimensiosn before removing cells")
print(dim(count.dat$counts))
count.dat$counts <- count.dat$counts[, cells.keep]
print("Dimensiosn before after removing cells")
print(dim(count.dat$counts))

prop.vec.filtered <- prop.vec[cells.keep]

# remove cells that have cellsizes less than counts_final
counts_new <- DropletUtils::downsampleMatrix(count.dat$counts, prop = prop.vec.filtered)

print(range(colSums(count.dat$counts)))
print(range(colSums(counts_new)))

count.dat$counts <- counts_new

save(count.dat, file = args$outfile)
print(jtime - Sys.time())
