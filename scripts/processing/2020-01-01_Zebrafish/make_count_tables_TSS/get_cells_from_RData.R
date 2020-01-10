# Jake Yeung
# get_cells_from_RData.R
# 2019-11-08
# DESCRIPTION
# 
#     Get good cells from .RData
# 
# FOR HELP
# 
#     Rscript get_cells_from_RData.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-11-08
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='.RData with count.dat$counts')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Good cells text file')
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

if (endsWith(args$infile, ".RData")){
    load(args$infile)  # count.dat$counts
} else {
    count.dat <- list()
    count.dat$counts <- readRDS(args$infile)
}

good.cells <- colnames(count.dat$counts)

# write to output
sink(args$outfile)
for (gc in good.cells){
  cat(gc, "\n")
}
sink()



