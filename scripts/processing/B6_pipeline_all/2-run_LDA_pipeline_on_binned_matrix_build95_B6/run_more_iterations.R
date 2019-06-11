# Jake Yeung
# 2c-run_more_iterations.R
# 2019-05-12
# DESCRIPTION
# 
#     More iterations 
# 
# FOR HELP
# 
#     Rscript 2c-run_more_iterations.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-05-12
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(topicmodels)
library(Matrix)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='LDA object output with multiple topics')
parser$add_argument('mark', metavar='MARK',
                                            help='Mark')
parser$add_argument('kchoose', metavar='K', type="integer", 
                                            help='Select topic')
parser$add_argument('niter', metavar='niter', type="integer", 
                                            help='Number of iterations')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='LDA object new')
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

jmark <- args$mark
kchoose <- args$kchoose

load(infile, v=T)

out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)

opts <- list(iter = args$niter)
lda.out <- LDA(count.mat, args$kchoose, method = "Gibbs", control = opts, model = out.objs$out.lda)

saveRDS(lda.out, file = args$outfile)
