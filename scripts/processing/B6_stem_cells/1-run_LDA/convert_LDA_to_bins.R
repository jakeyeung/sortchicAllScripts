# Jake Yeung
# 3-convert_LDA_to_bins.R
# 2019-09-30
# DESCRIPTION
# 
#     Description
# 
# FOR HELP
# 
#     Rscript 3-convert_LDA_to_bins.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-09-30
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(scchicFuncs)
library(topicmodels)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input LDA Output')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output LDA bins')
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

load(args$infile, v=T)  # out.lda

out.objs <- LoadLDABins(jmark = "H3K4me1", inf = args$infile, choose.k = 50)

save(out.objs, file = args$outfile)
