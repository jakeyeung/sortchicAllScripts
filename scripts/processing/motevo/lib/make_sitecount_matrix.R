# Jake Yeung
# 5-make_sitecount_matrix.R
# 2019-02-04
# DESCRIPTION
# 
#     Make sitecount matrix from sql db
# 
# FOR HELP
# 
#     Rscript 5-make_sitecount_matrix.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-02-04
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time()

library(dplyr)
library(data.table)
library(tidyr)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='SQL DB file')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Sitecount matrix ready for MARA')
parser$add_argument('--wmlist', metavar='WMs List', required = TRUE,
                                            help='Get list of WMs. Path to file')
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

WMs <- unlist(data.table::fread(args$wmlist, header=FALSE))
# for version 2, take out .pwm
# WMs <- gsub(".pwm", "", WMs)

names(WMs) <- WMs

motevo.db <- src_sqlite(args$infile, create=F)

tblname <- src_tbls(motevo.db)[[1]]

motevo.tbl <- tbl(motevo.db, sql(paste0("SELECT * FROM ", tblname)))

dat.lst <- lapply(WMs, function(x) motevo.tbl %>% filter(motif == x) %>% group_by(peak, motif) %>% summarise(nsites = length(sitecount), sitecount = sum(sitecount)))

dat.raw <- lapply(dat.lst, function(x) collect(x, n = Inf))

dat.long <- dplyr::bind_rows(dat.raw)

dat.mat <- as.data.frame(tidyr::spread(subset(dat.long, select = -c(nsites)), motif, sitecount))
Gene.ID <- dat.mat$peak
rownames(dat.mat) <- dat.mat$peak; dat.mat$peak <- NULL

dat.mat <- as.matrix(dat.mat)
dat.mat[is.na(dat.mat)] <- 0

# normalize
dat.mat <- t(scale(t(scale(dat.mat, center=TRUE, scale=TRUE)), center=TRUE, scale=TRUE))

dat.mat <- as.data.frame(dat.mat)
dat.mat <- cbind(Gene.ID, dat.mat)
# rownames(dat.mat) <- jpeaks

print(head(dat.mat))
# write to output
data.table::fwrite(dat.mat, file = args$outfile, sep="\t")

print(Sys.time() - jstart)
