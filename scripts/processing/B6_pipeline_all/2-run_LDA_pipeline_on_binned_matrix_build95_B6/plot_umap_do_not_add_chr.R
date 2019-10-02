# Jake Yeung
# plot_umap.R
# 2019-05-08
# DESCRIPTION
# 
#     Load LDA object, plot UMAP
# 
# FOR HELP
# 
#     Rscript plot_umap.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-05-08
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")

tstart <- Sys.time()

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Input .Rdata object of LDA')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output pdf')
parser$add_argument('outrdata', metavar='OUTDATA',
                                            help='Output rdata')
parser$add_argument('mark', metavar='MARK',
                                            help='Mark: H3K4me1 for example')
parser$add_argument('thres', metavar='THRES', type="double", 
                                            help='Threshold for significant peak')
parser$add_argument('ncores', metavar='NCORES', type="integer", 
                                            help='Number of cores')
parser$add_argument('nnmin', metavar='nnmin', type="integer", 
                                            help='nn min')
parser$add_argument('nnmax', metavar='nnmax', type="integer", 
                                            help='nn max')
parser$add_argument('nnby', metavar='nnby', type="integer", 
                                            help='nn by')
parser$add_argument('mindist', metavar='mindsit', type="double", 
                                            help='mindist umap parameter')
parser$add_argument('kchoose', metavar='kchoose', type="character", 
                                            help='k choose')
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

jseed.louv=123
jmetric.louv='euclidean'
jmindist <- args$mindist
jnn.vec <- seq(args$nnmin, args$nnmax, by = args$nnby)
names(jnn.vec) <- as.character(jnn.vec)
print("Iterating through:")
print(jnn.vec)

# custom.settings.new.lst <- lapply(jnn.vec, function(x) GetUmapSettings(x, jmetric.louv, jmindist, jseed.louv))

print("Loading bins...")
out.objs <- LoadLDABins(args$mark, jbin = NA, top.thres = args$thres, inf = args$infile, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = FALSE, choose.k = args$kchoose)
kchoose.str <- paste0("K:", args$kchoose, "_", out.objs$out.lda@k)
print(kchoose.str)

print("Making umaps")
# for (nn in jnn.vec){
mlst <- parallel::mclapply(jnn.vec, function(nn){
# for (nn in jnn.vec){
  custom.settings <- GetUmapSettings(nn, jmetric.louv, jmindist, jseed.louv)
  dat.umap <- umap(out.objs$tm.result$topics, config = custom.settings)
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2])
  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste("NN:", nn, "metric:", jmetric.louv, "MinDist:", jmindist, "K:", kchoose.str))
  return(m)
# }
}, mc.cores = args$ncores)

print("Plotting to output")
pdf(args$outfile, useDingbats=FALSE)
  for (nn in jnn.vec){
    print(paste("Plotting topic:", nn))
    print(mlst[[as.character(nn)]])
  }
dev.off()
# save Rdata to ouptput
save(out.objs, mlst, jnn.vec, jmetric.louv, jmindist, jseed.louv, file = args$outrdata)
print(Sys.time() - tstart)
