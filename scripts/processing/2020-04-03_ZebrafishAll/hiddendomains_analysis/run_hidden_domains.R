# Jake Yeung
# run_hidden_domains.R
# 2020-04-24
# DESCRIPTION
# 
#     Run hidden domains
#     Trouble shooting tips on outputs with no HDs called on a chromosmoe:  
#     http://hiddendomains.sourceforge.net
# 
# FOR HELP
# 
#     Rscript run_hidden_domains.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-24
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
source("/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/hiddenDomains.R")  # hiddenDomains funcion

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE',
                                            help='Input treatments_bins')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output file .txt')
parser$add_argument('-maxreadcount', metavar='Whole number', type="integer", default=10,
                                            help='Max read count, set to 10 if HD misses chromosomes. Otherwise 200?')
parser$add_argument('-minreadcount', metavar='Whole number', type="integer", default=-1,
                                            help='Can fiddle if HMM fails to conerge, set to -1 if things dont work otherwise -10?')
parser$add_argument('-minprob', metavar='Prob', type="double", default=0.6,
                                            help='Filters out HD less than minprob')
parser$add_argument('-chromonames', metavar='space separated strings', nargs="+", default="mouse",
                                            help='Insert vector of chromonames. mouse or human also accepted')
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

# returns nothing?
hiddenDomains(treat.bin.file = args$infile, control.bin.file=NULL,  min.prob=args$minprob, max.read.count = args$maxreadcount, min.read.count = args$minreadcount, out.file.name = args$outfile, chr.names=args$chromonames)

