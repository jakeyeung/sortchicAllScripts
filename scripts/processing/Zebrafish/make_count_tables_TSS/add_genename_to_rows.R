# Jake Yeung
# add_genename_to_rows.R
# 2019-11-11
# DESCRIPTION
# 
#     Add genename to rows
# 
# FOR HELP
# 
#     Rscript add_genename_to_rows.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-11-11
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(data.table)
library(hash)
library(Matrix)

GetGeneAnnotsHash <- function(inf.annot){
  dat.annot <- data.table::fread(inf.annot, col.names = c("chromo", "start", "end", "bname"))
  # add chr
  dat.annot$chromo <- paste("chr", dat.annot$chromo, sep = "")
  rnames.old <- paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":")
  rnames.new <- dat.annot$bname
  annots.hash <- hash::hash(rnames.old, rnames.new)
}

AddGeneNameToRows <- function(mat, annots.hash){
  # mat rownmaes got stripped of gene names, add them back
  rnames.old <- rownames(mat)
  rnames.new <- sapply(rnames.old, function(x) annots.hash[[x]])
  rownames(mat) <- rnames.new
  return(mat)
}

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='.RData containing count.dat$counts')
parser$add_argument('infannot', metavar='ANNOT',
                                            help='Bed file with 4th column containing new rownames')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='.RData containing count.dat$counts new rownames')
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

if (endsWith(args$infile, ".rds")){
  count.dat <- list()
  count.dat$counts <- readRDS(args$infile)
} else {
  load(args$infile, v=T)  # count.dat$counts
}

annots.hash <- GetGeneAnnotsHash(args$infannot)

rnames.new <- sapply(rownames(count.dat$counts), function(x){
  xnew <- annots.hash[[x]]                       
  if (is.null(xnew)){
    print(paste("Warning: x does not match anything in hash, returning orig: ", x))
    return(x)
  } else {
    return(xnew)
  }
})

rownames(count.dat$counts) <- rnames.new

print("Names before: ")
print(head(rnames.new))
print("Names after: ")
print(head(rownames(count.dat$counts)))

# write output
save(count.dat, file = args$outfile)
