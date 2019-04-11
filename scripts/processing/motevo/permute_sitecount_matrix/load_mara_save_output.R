# Jake Yeung
# load_mara_save_output.R
# 2019-04-10
# DESCRIPTION
# 
#     Load MARA using my script then save output
# 
# FOR HELP
# 
#     Rscript load_mara_save_output.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-10
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(hash)
library(data.table)
library(dplyr)
library(tidyr)
# setwd("/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic")
setwd("/home/hub_oudenaarden/jyeung/projects/scChiC")
source("scripts/Rfunctions/MatchCellNameToSample.R")
source("scripts/Rfunctions/MaraDownstream.R")

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('indir', metavar='INDIR',
                                            help='indir')
parser$add_argument('mark', metavar='MARK',
                                            help='Mark for labeling column name')
parser$add_argument('seed', metavar='SEED',
                                            help='Seed for annotating the permutation seed')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='outfile')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")

                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

assertthat::assert_that(!file.exists(args$outfile))

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

mdir <- args$indir

experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
switch.rep.hash <- GetRepSwitchHash(experihash)

mara.out <- LoadMARA(mdir, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = ".gz")
mara.out$act.long$mark <- args$mark
mara.out$act.long$seed <- args$seed
data.table::fwrite(mara.out$act.long, file = args$outfile, sep = "\t")

system(paste("gzip", args$outfile))
