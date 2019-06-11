# Jake Yeung
# merge_count_matrix.R
# 2019-04-19
# DESCRIPTION
# 
#     Merge two count mats from Robj
# 
# FOR HELP
# 
#     Rscript merge_count_matrix.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-19
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile1', metavar='RData path',
                                            help='Rdata 1')
parser$add_argument('infile2', metavar='RData path',
                                            help='Rdata 2')
parser$add_argument('outfile', metavar='Rdata path',
                                            help='Rdata output')
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

load(args$infile1, v=T)
mat1 <- count.dat$counts 
load(args$infile2, v=T)
mat2 <- count.dat$counts

mat.merge <- rbind(mat1, mat2)

count.dat$counts <- mat.merge

save(count.dat, file = args$outfile)
