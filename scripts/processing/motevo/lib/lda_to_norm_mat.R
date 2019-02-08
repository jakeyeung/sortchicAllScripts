# Jake Yeung
# lda_to_norm_mat.R
# 2019-02-04
# DESCRIPTION
# 
#     Load LDA output and write normalized count matrix 
# 
# FOR HELP
# 
#     Rscript lda_to_norm_mat.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-02-04
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

library(dplyr)
library(topicmodels)
library(data.table)

suppressPackageStartupMessages(library("argparse"))

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='.RData object containing out.lda')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Text output')
parser$add_argument("-t", "--thres", type = "double", default=0.995,
                        help="Threshold for picking top hits")
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

load(args$infile, v=T)  # out.lda, counts.mat

if (length(out.lda) > 1){
  out.lda <- ChooseBestLDA(out.lda)
}
kchoose <- out.lda@k

tm.result <- topicmodels::posterior(out.lda)

topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = args$thres))
})
top.regions <- unique(unlist(topic.regions))
mat.norm <- as.data.frame(t(tm.result$topic %*% tm.result$terms))
print(dim(mat.norm))
mat.norm <- mat.norm[top.regions, ]
print(dim(mat.norm))


# log? 
mat.norm <- log2(mat.norm * 10^6 + 1)

Gene.ID <- rownames(mat.norm)

dat.out <- cbind(Gene.ID, mat.norm)

data.table::fwrite(dat.out, file=args$outfile, sep="\t", row.names=FALSE)

