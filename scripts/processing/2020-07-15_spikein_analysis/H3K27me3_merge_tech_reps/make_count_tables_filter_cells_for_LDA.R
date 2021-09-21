# Jake Yeung
# make_count_tables_filter_cells_for_LDA.R
# 2020-11-28
# DESCRIPTION
# 
#     Load .csv count data, filter out cells based on annot data (column cell)
# 
# FOR HELP
# 
#     Rscript make_count_tables_filter_cells_for_LDA.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-11-28
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

cbind.fill.lst <- function(mats.lst, all.rnames, fill = 0){
  mats.lst.filled <- lapply(mats.lst, function(mat.tmp){
    missing.rnames <- all.rnames[!all.rnames %in% rownames(mat.tmp)]
    mat.tmp.to.fill <- matrix(data = fill, nrow = length(missing.rnames), ncol = ncol(mat.tmp), dimnames = list(missing.rnames, colnames(mat.tmp)))
    mat.tmp.bind <- rbind(mat.tmp, mat.tmp.to.fill)
    mat.tmp.bind <- mat.tmp.bind[all.rnames, ]
    return(mat.tmp.bind)
  })
  return(do.call(cbind, mats.lst.filled))
}


suppressPackageStartupMessages(library("argparse"))
library(scchicFuncs)
library(Matrix)
library(dplyr)
library(data.table)
library(assertthat)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE', nargs='+',
                                            help='Input csv file, can be a vector of more than one file then willi be cbinded with fill')
parser$add_argument('-annotfile', metavar='INFILE',
                                            help='Input csv file, with colum name "cell" to specify good cells')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output .rds file for LDA')
parser$add_argument("--add_chromo", action="store_true",
                        help="Add chr prefix to rownames")
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

dat.meta <- fread(args$annotfile)
good.cells <- dat.meta$cell
assertthat::assert_that(length(good.cells) > 0)

print(paste("Number of input files:", length(args$infile)))
if (length(args$infile) == 1){
  print("Reading single file")
  mat <- ReadMatSlideWinFormat(inf = args$infile, as.sparse = TRUE, add.chromo = args$add_chromo)
} else {
  print("Reading multiple files and cbinding with fill")
  infs.lst <- as.list(args$infile)
  names(infs.lst) <- args$infile
  mat.lst <- lapply(infs.lst, function(inftmp){
    mat <- ReadMatSlideWinFormat(inf = inftmp, as.sparse = TRUE, add.chromo = args$add_chromo)
  })
  all.rnames <- unique(unlist(lapply(mat.lst, rownames)))
  rnames.i <- gtools::mixedorder(all.rnames)
  all.rnames.sorted <- all.rnames[rnames.i]
  mat <- cbind.fill.lst(mat.lst, all.rnames = all.rnames.sorted)
}
print("Mat loaded of size:")
print(dim(mat))

colskeep <- which(colnames(mat) %in% good.cells)

print("Filtering good cells:")
print(length(good.cells))
assertthat::assert_that(length(colskeep) > 0)

mat.filt <- mat[, colskeep]

saveRDS(mat.filt, file = args$outfile)


