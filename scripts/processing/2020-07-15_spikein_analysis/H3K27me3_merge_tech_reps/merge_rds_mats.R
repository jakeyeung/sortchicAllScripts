# Jake Yeung
# merge_rds_mats.R
# 2020-11-28
# DESCRIPTION
# 
#     Take multipe rds objects and cbind them 
# 
# FOR HELP
# 
#     Rscript merge_rds_mats.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-11-28
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(Matrix)

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

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE', nargs="+",
                                            help='List of rds object to cbind')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output .rds file ready for LDA')
# parser$add_argument("--common_rows", action="store_true",
#                         help="Take common rows, otherwise take union")
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

print("Reading multiple files and cbinding with fill")
infs.lst <- as.list(args$infile)
names(infs.lst) <- args$infile
mat.lst <- lapply(infs.lst, function(inftmp){
  print(inftmp)
  mattmp <- readRDS(inftmp)
  print(dim(mattmp))
  return(mattmp)
})
all.rnames <- unique(unlist(lapply(mat.lst, rownames)))
common.rnames <- Reduce(f = intersect, x = lapply(mat.lst, rownames))
print(paste("Union rownames:", length(all.rnames)))
print(paste("INtersect rownames:", length(common.rnames)))
rnames.i <- gtools::mixedorder(all.rnames)
all.rnames.sorted <- all.rnames[rnames.i]
mat <- cbind.fill.lst(mat.lst, all.rnames = all.rnames.sorted)

saveRDS(mat, file = args$outfile)

