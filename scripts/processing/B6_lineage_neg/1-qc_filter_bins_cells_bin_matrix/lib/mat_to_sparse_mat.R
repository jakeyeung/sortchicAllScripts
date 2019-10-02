# Jake Yeung
# mat_to_sparse_mat.R
# 2019-05-07
# DESCRIPTION
# 
#     Mat to Sparse Mat
# 
# FOR HELP
# 
#     Rscript mat_to_sparse_mat.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-05-07
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(Matrix)
library(data.table)
library(dplyr)

ReadMat <- function(inf, as.sparse = TRUE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  startend = V2)
  # remove missing
  dat <- subset(dat, startend != "Missing")
  dat$coord <- paste(dat$chromo, dat$startend, sep = ":")
  dat$chromo <- NULL
  dat$startend <- NULL
  if (as.sparse){
    coords <- dat$coord
    dat$coord <- NULL
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
    rownames(dat) <- coords
    dat <- Matrix::Matrix(dat, sparse = TRUE)
  }
  return(dat)
}

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Infile')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Out RDS')
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

mat <- ReadMat(args$infile, as.sparse = TRUE)
saveRDS(mat, file = args$outfile)


