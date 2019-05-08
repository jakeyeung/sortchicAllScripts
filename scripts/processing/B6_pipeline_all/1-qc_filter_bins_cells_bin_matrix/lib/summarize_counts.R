# Jake Yeung
# summarize_counts.R
# 2019-05-07
# DESCRIPTION
# 
#     Summarize counts
# 
# FOR HELP
# 
#     Rscript summarize_counts.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-05-07
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(dplyr)
library(tidyr)
library(data.table)

GetMarkB6 <- function(x){
  return(strsplit(x, "-")[[1]][[4]])
}

GetRepB6 <- function(x){
  return(paste0("rep", strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
}

GetCellB6 <- function(x, shift = 1){
  # shift by 1 to get cells to start at 1 rather than 0
  cell <- as.numeric(strsplit(x, "_")[[1]][[4]]) + shift
  return(as.character(cell))
}

ReadDinuc <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName) 
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc) %>%
    rowwise() %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkB6(samp),
           repl = GetRepB6(samp),
           cell = GetCellB6(samp))
  return(dat.long)
}

ReadMat <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  startend = V2)
  # remove missing
  dat <- subset(dat, startend != "Missing")
  dat$coord <- paste(dat$chromo, dat$startend, sep = ":")
  dat$chromo <- NULL
  dat$startend <- NULL
  return(dat)
}


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input file')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Summaized across datasets')
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

dat.mat <- ReadMat(args$infile)
dat.cellsums <- colSums(dat.mat %>% dplyr::select(-coord), na.rm = TRUE)
dat.out <- data.frame(cell = names(dat.cellsums), cellsum = dat.cellsums)

# write output
fwrite(dat.out, file = args$outfile)

# system(paste0("gzip ", args$outfile))
