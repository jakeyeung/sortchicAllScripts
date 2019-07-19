# Jake Yeung
# get_best_LDA_from_list.R
# 2019-06-18
# DESCRIPTION
# 
#     Get best LDA
# 
# FOR HELP
# 
#     Rscript get_best_LDA_from_list.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-18
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input RData of list of LDA objects')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output RDS of best LDA object, or pick k as you want')
parser$add_argument("--best_k", default="auto", help="Manually choose K from list of LDA")
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


