# Jake Yeung
# 1-merge_count_tables_into_one.R
# 2021-02-15
# DESCRIPTION
# 
#     Merge count tables into one
# 
# FOR HELP
# 
#     Rscript 1-merge_count_tables_into_one.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2021-02-15
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

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

GetCountMat <- function(infs.k4){
  # (infs.k4 <- list.files(indir, pattern = "BM_round1_round2_merged_H3K4me1_.*..bam.count_table_k4_k9_dynamic_regions.txt", full.names = TRUE))
  counts.mat.k4.lst <- lapply(infs.k4, function(inf){
    count.mat.tmp <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
    return(count.mat.tmp)
  })
  rnames.all.k4 <- unique(unlist(lapply(counts.mat.k4.lst, function(x) rownames(x))))
  count.mat.k4 <- cbind.fill.lst(counts.mat.k4.lst, all.rnames = rnames.all.k4, fill = 0)
}


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE', nargs="+", 
                                            help='Space delimited list of input mats')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output rds')
parser$add_argument('-metafile', metavar='METAFILE',
                                            help='metafile')
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


# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# # hubprefix <- "/home/jyeung/hub_oudenaarden"
# hubprefix <- "/hpc/hub_oudenaarden"

dat.meta <- fread(args$metafile)

cells.keep <- dat.meta$cell

# Get singles -------------------------------------------------------------

infs.lst <- args$infile
count.mat <- GetCountMat(infs.lst)
print(dim(count.mat))
print("remove duplicate columns?")
count.mat <- count.mat[, !duplicated(colnames(count.mat))]
print(dim(count.mat))
print("Keep good cells")
cells.keep.i <- colnames(count.mat) %in% cells.keep
count.mat <- count.mat[, cells.keep.i]
print(dim(count.mat))
saveRDS(count.mat, file = args$outfile)

