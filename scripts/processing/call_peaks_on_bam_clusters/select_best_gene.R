# Jake Yeung
# select_best_gene.R
# 2019-04-19
# DESCRIPTION
# 
#     Select best gene
# 
# FOR HELP
# 
#     Rscript select_best_gene.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-04-19
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(dplyr)
library(tidyr)
library(Matrix)

tstart <- Sys.time()


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Select best gene out of many isoforms based on total reads across cells')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='New Robj')
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

load(args$infile, v=T)

# rectangularize because I can

gvec <- sapply(rownames(count.dat$counts), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
print("my gvec")
print(gvec[1:5])
dat.long <- data.frame(pname = rownames(count.dat$counts), gene = gvec, as.data.frame(as.matrix(count.dat$counts))) %>% gather(key = "cell", value = "counts", c(-pname, -gene))

print(head(dat.long))

print("Summarizing...")
# now we can select
dat.sum <- dat.long %>% 
  group_by(pname, gene) %>% 
  summarise(counts = sum(counts)) %>%
  group_by(gene) %>% 
  filter(counts == max(counts)) %>% 
  filter(counts > 0) %>% 
  mutate(ngenes = length(counts))
print("Summarizing...done")

print(head(dat.sum))
dat.sum <- dat.sum[!duplicated(dat.sum$gene), ]

# regenerate the count matrix with new filtered rownames

peaks.filt <- dat.sum$pname

row.i <- rownames(count.dat$counts) %in% peaks.filt
cmat.new <- count.dat$counts[row.i, ]

count.dat$counts <- cmat.new

save(count.dat, dat.sum, file = args$outfile)

print(Sys.time() - tstart)
