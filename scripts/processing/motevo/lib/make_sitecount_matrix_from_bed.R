# Jake Yeung
# make_sitecount_matrix_from_bed.R
# 2019-02-04
# DESCRIPTION
# 
#     Make sitecount matrix from bed. Loads up 20 GBs of data??
# 
# FOR HELP
# 
#     Rscript make_sitecount_matrix_from_bed.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-02-04
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time()

library(data.table)
library(dplyr)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='In bed ~20 GB? 8 Columns. No strand')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Out sitecount matrix')
parser$add_argument('--scale', metavar="TRUE or FALSE", type = "integer", 
                    default=0, help="Scale matrix?")
parser$add_argument('--center', metavar="TRUE or FALSE", type = "integer", 
                    default=0, help="Center matrix?")
parser$add_argument('--byrow', action="store_true", 
                    default=FALSE, help="Also normalize matrix by row")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# change int to logical
args$scale <- as.logical(args$scale)
args$center <- as.logical(args$center)

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

dat <- fread(args$infile, header=FALSE)
if (ncol(dat) == 8){
  colnames(dat) <- c("chromo", "start", "end", "motif", "sitecount", "gene", "dist", "peak")
}

dat.sum <- dat %>% group_by(peak, motif) %>% summarise(nsites = length(sitecount), sitecount = sum(sitecount))

dat.mat <- as.data.frame(tidyr::spread(subset(dat.sum, select = -c(nsites)), motif, sitecount, fill=0))

Gene.ID <- dat.mat$peak
rownames(dat.mat) <- dat.mat$peak; dat.mat$peak <- NULL

# normalize?
# dat.mat <- t(scale(t(scale(dat.mat, center=TRUE, scale=TRUE)), center=TRUE, scale=TRUE))
dat.mat <- scale(dat.mat, center=args$center, scale=args$scale)

if (args$byrow){
  dat.mat <- t(scale(t(dat.mat), center=args$center, scale=args$scale))
}

dat.mat <- as.data.frame(dat.mat)
dat.mat <- cbind(Gene.ID, dat.mat)

print(head(dat.mat))
# write to output
data.table::fwrite(dat.mat, file = args$outfile, sep="\t")

print(Sys.time() - jstart)

