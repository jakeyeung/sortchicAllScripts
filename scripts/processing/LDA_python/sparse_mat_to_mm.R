# Jake Yeung
# sparse_mat_to_mm.R
# 2019-07-11
# DESCRIPTION
# 
#     Read sparse mat as RData, write to mm (expect count.dat$counts)
# 
# FOR HELP
# 
#     Rscript sparse_mat_to_mm.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-07-11
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(Matrix)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input RData object with count.dat$counts')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='File output to mm')
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

writeMM(count.dat$counts, sparse = TRUE, file = args$outfile)

# write rownames and column names
write.table(rownames(count.dat$counts), file = paste0(args$outfile, ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(colnames(count.dat$counts), file = paste0(args$outfile, ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
