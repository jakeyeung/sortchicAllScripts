# Jake Yeung
# filter_correlated_bins_across_mats.R
# 2020-04-08
# DESCRIPTION
# 
#     remove correlated plots across multiple marks 
# 
# FOR HELP
# 
#     Rscript filter_correlated_bins_across_mats.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-04-08
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(assertthat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(scchicFuncs)
library(Matrix)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(tidyr)
library(GGally)
library(here)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infiles', metavar='INFILE', nargs='+',
                                            help='List of input .rds files across multiple marks')
parser$add_argument('-outfiles', metavar='OUTFILE', nargs='+',
                                            help='Output .rds files, filtering out correlated bins. Should match 1to1 with infiles')
parser$add_argument('-names', metavar='HistoneMod', nargs='+',
                                            help='Names associated with each .rds file. Should match 1to1 with infiles')
parser$add_argument('-quantilethres', metavar='Quantile', type='double', default=0.95,
                                            help='Quantile threshold to remove correlated bins')
parser$add_argument('-badbinsoutbed', metavar='BEDtxt', 
                                            help='Output bed file of bad correlated bins that were removed')
parser$add_argument('-pdfout', metavar='OUTF', 
                                            help='Output pdf of bins that were removed')
parser$add_argument('--overwrite', default=FALSE, action='store_true',
                                            help='Overwrites existing files ')
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

infiles <- args$infiles
outfiles <- args$outfiles
qthres <- args$quantilethres
jmarks <- args$names

badbins.outbed <- args$badbinsoutbed
pdfout <- args$pdfout

if (!args$overwrite){
  assertthat::assert_that(!file.exists(badbins.outbed))
  assertthat::assert_that(!file.exists(pdfout))
  lapply(outfiles, function(outf){
    assertthat::assert_that(!file.exists(outf))
  })
}

names(infiles) <- jmarks
names(outfiles) <- jmarks

assertthat::assert_that(length(infiles) == length(outfiles))

mats.lst <- lapply(infiles, function(inf){
  mat <- readRDS(inf)
  return(mat)
})

bin.lst <- lapply(mats.lst, function(mat){
  b <- Matrix::rowMeans(mat)
})

rnames.lst <- lapply(bin.lst, function(b){
  names(b)
})

print(paste("Number of bins for each mark:"))
lengths.all <- lapply(rnames.lst, function(x) length(x))
print(lengths.all)

rnames.all.common <- Reduce(intersect, rnames.lst)  # intersect
print(paste("Number of common bins before filtering:"))
print(length(rnames.all.common))

bin.long <- lapply(jmarks, function(jmark){
  x <- bin.lst[[jmark]]
  dat <- data.frame(coord = names(x), bin.mean = x, mark = jmark, stringsAsFactors = FALSE)
  dat <- subset(dat, coord %in% rnames.all.common)
}) %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(bin.rank = rank(bin.mean) / length(bin.mean))  # ranks are normalized 0 to 1
bin.wide <- spread(bin.long %>% dplyr::select(-bin.rank), key = mark, value = bin.mean)

bin.sum <- bin.long %>%
  group_by(coord) %>%
  summarise(bin.rank.mean = mean(bin.rank),
            bin.rank.min = min(bin.rank),
            bin.rank.max = max(bin.rank))

# remove bad bins?
bins.high <- subset(bin.sum, bin.rank.mean >= qthres)$coord
bins.low <- subset(bin.sum, bin.rank.mean <= (1-qthres))$coord

bins.correlated <- c(bins.high, bins.low)
print("Number of correlated bins to remove:")
print(length(bins.correlated))

# plot output

# m.pairs.clean <- ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
#                          lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
#   ggtitle("After removing correlated bins")

# faster to filter on the wide bins
print("Saving correlated bins to variable")
print(head(bin.wide))
rnames.corrfilt <- (bin.wide %>% dplyr::filter(!coord %in% bins.correlated))$coord

print("Filtering bins for each mat")
mats.filt.lst <- lapply(mats.lst, function(mat){
  rnames.keep <- rownames(mat) %in% rnames.corrfilt
  mat.filt <- mat[rnames.keep, ]
})

print("Mat dimensions before:")
lapply(mats.lst, dim)

print("Mat dimensions after:")
lapply(mats.filt.lst, dim)

# write to output
lapply(jmarks, function(jmark){  print(jmark)
  mat <- mats.filt.lst[[jmark]]
  print(jmark)
  print(dim(mat))
  print(1 - Matrix::nnzero(mat) / length(mat))
  saveRDS(mat, file = outfiles[[jmark]])
})

# write correlated bins to output
bins.correlated.dat <- data.frame(chromo = paste0("chr", sapply(bins.correlated, GetChromo)),
                                  start = sapply(bins.correlated, GetStart),
                                  end = sapply(bins.correlated, GetEnd))
fwrite(bins.correlated.dat, file = badbins.outbed, col.names = FALSE, sep = "\t")

# write PDFs 
pdf(file = pdfout, useDingbats=FALSE)
  ggpairs(log10(bin.wide %>% dplyr::filter(!coord %in% bins.correlated) %>% dplyr::select(-coord)),
          lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
    ggtitle("After removing correlated bins")
  
  ggpairs(log10(bin.wide %>% dplyr::select(-coord)), lower = list(continuous = wrap("points", alpha = 0.1))) + theme_classic() +
    ggtitle("Before removing correlated bins")
dev.off()

