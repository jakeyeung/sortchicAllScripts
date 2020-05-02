# Jake Yeung
# plot_density_treatment_bins.R
# 2020-04-24
# DESCRIPTION
# 
#     Find maxcounts threshold for calling hiddenDomains
# 
# FOR HELP
# 
#     Rscript plot_density_treatment_bins.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-24
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(data.table)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input treatment_bins.txt')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output .pdf')
parser$add_argument('-threshold', metavar='Number', type="double", default=0,
                                            help='Draw vertical threshold')
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

dat <- fread(args$infile)  # colnames: id      chr     pos     count

pdf(args$outfile, useDingbats=FALSE)
plot(density(dat$count))
abline(v = args$threshold)
hist(dat$count)
abline(v = args$threshold)
plot(density(dat$count), log="x")
abline(v = args$threshold)
hist(log10(dat$count + 1))
abline(v = log10(args$threshold + 1))
dev.off()
