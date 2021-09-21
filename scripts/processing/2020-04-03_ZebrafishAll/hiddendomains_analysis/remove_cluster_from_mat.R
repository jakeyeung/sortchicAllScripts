# Jake Yeung
# remove_cluster_from_mat.R
# 2020-08-25
# DESCRIPTION
# 
#     Remove cells assigned to cluster
# 
# FOR HELP
# 
#     Rscript remove_cluster_from_mat.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-08-25
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(data.table)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE',
                                            help='Input .rds')
parser$add_argument('-annotfile', metavar='INFILE',
                                            help='Input .txt with cluster column name')
parser$add_argument('-clstremove', metavar='Celltype',
                                            help='Name of cluster to remove')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output .rds with cells removed')
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

mat <- readRDS(args$infile)

dats.annot <- fread(args$annotfile)

# cells.keep <- subset(dats.annot, cluster != args$clstremove)$cell
cells.keep <- subset(dats.annot, !grepl(args$clstremove, cluster))$cell

cells.keep.i <- colnames(mat) %in% cells.keep

assertthat::assert_that(length(cells.keep) > 0)

mat.filt <- mat[, cells.keep.i]

saveRDS(mat.filt, file = args$outfile)

