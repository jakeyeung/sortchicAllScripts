# Jake Yeung
# 1-process_zf_scrnaseq.R
# 2019-08-16
# DESCRIPTION
# 
#     Read csv.gz scrnaseq zebrafisha nd output, filter cells make rds ready for transform
# 
# FOR HELP
# 
#     Rscript 1-process_zf_scrnaseq.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-08-16
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(Matrix)
library(data.table)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input csv.gz file the_massive_complete_zf_dataset.csv.gz. Columns are cells, rows are genes. No column for rownames')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output rds matrix')
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

dat <- fread(args$infile)
colnames(dat)[[1]] <- "gene"

mat <- Matrix(as.matrix(as.data.frame(dat[, -1])), sparse = TRUE)
rownames(mat) <- dat$gene

cells.keep <- which(colSums(mat) > 1000)

saveRDS(mat[, cells.keep], file = args$outfile)



