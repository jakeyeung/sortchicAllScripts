# Jake Yeung
# project_new_samples_on_LDA.R
# 2019-06-18
# DESCRIPTION
# 
#     Project new samples (in form of a sparse matrix) onto existing LDA
# 
# FOR HELP
# 
#     Rscript project_new_samples_on_LDA.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-18
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(topicmodels)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE RDATA',
                                            help='Input .RData containing LDA object')
parser$add_argument('inmat', metavar='INMAT RDS',
                                            help='Input .rds containing sparse matrix (bins in rows, samples in columsn)')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output RData containing posterior of predicted counts')
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


load(args$infile, v=T)  # out.objs$out.lda 
assertthat::assert_that(class(out.objs$out.lda) == "LDA_Gibbs")

load(args$inmat, v=T)  # count.dat$counts
assertthat::assert_that(class(count.dat$counts) == "dgCMatrix")

system.time(
  out.lda.predict <- posterior(out.objs$out.lda, t(as.matrix(count.dat$counts)))
)

save(out.lda.predict, file = args$outfile)

