# Jake Yeung
# run_sctransform.R
# 2019-08-16
# DESCRIPTION
# 
#     Run sctransform
# 
# FOR HELP
# 
#     Rscript run_sctransform.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-08-16
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(Matrix)
library(Seurat)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input .rds file')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output Seurat object rds file')
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

mat <- readRDS(args$infile)  # rows are genes, columns are cells, filter out cells beforehand

sobj <- CreateSeuratObject(counts = mat)
sobj <- SCTransform(sobj, verbose=TRUE)

saveRDS(sobj, file = args$outfile)

