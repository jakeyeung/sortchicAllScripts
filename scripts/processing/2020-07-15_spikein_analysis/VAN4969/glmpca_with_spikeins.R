# Jake Yeung
# glmpca_with_spikeins.R
# 2020-08-11
# DESCRIPTION
# 
#     Run glmpca with spikeins

# 
# FOR HELP
# 
#     Rscript glmpca_with_spikeins.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2020-08-11
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 


library(scchicFuncs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument('-infile', metavar='INFILE',
                    help='Input .rds file count matrix')
parser$add_argument('-outfile', metavar='OUTFILE',
                    help='Output glmpca RData')
parser$add_argument('-inspike', metavar='Input',
                    help='RData containing dat.spikeins.mat object')
parser$add_argument('-K', metavar='Dimensions', type = 'integer', default = 30,
                    help='Number of dimensions')
parser$add_argument('-penalty', metavar='Penalization term', type = 'double', default = 1.0,
                    help='Penalty (numeric)')
parser$add_argument('-maxIter', metavar='maxIter', type = 'integer', default = 1000,
                    help='Max iterations. Default 1000')
parser$add_argument('-minibatch', metavar='minibatch', default = 'none',
                    help='Minibatch: c("none", "stochastic", "memoized") default none')
parser$add_argument('-optimizer', metavar='optimizer', default = 'fisher',
                    help='Optimizer: c("avagrad", "fisher") default none')
parser$add_argument('-tol', metavar='tolerance', default = 1e-4, type = 'double',
                    help='Tolerance default 1e-4')
parser$add_argument('-topndevgenes', metavar='TopNGenesByDeviance', default = 0, type = 'integer',
                    help='Filter genes by top, default 0 takes all genes')
parser$add_argument("--by_plate", action="store_true",
                    help="Add design matrix by platement")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

args <- parser$parse_args()
assertthat::assert_that(file.exists(args$infile))
assertthat::assert_that(file.exists(args$inspike))
assertthat::assert_that(!file.exists(args$outfile))


# Load spikeincounts ------------------------------------------------------

# 
# hubprefix <- "/home/jyeung/hub_oudenaarden"
# indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
# infs.chromo <- list.files(indir, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
# 
# jspikeinchromo <- "J02459.1"
# jchromos <- paste("", seq(19), sep = "")
# 
# dat.chromos <- lapply(infs.chromo, function(inf){
#   dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = NA)
# }) %>%
#   bind_rows() %>%
#   filter(chromo == "1")
# 
# dat.spikeins.mat <- as.data.frame(subset(dat.chromos, select = c(samp, spikeincounts)))
# rownames(dat.spikeins.mat) <- dat.spikeins.mat$samp
# 
# inspike <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData"
# save(dat.spikeins.mat, file = inspike)

if (endsWith(x = args$inspike, suffix = ".RData")){
  load(args$inspike, verbose=TRUE)  # dat.spikeins.mat
} else if (endsWith(x = args$inspike, suffix = ".rds")){
  dat.spikeins.mat <- readRDS(args$inspike)
  print("Loading from rds rownames shoudl be samp names")
  print(head(rownames(dat.spikeins.mat)))
} else if (endsWith(x = args$inspike, suffix = ".txt")){
  dat.spikeins.mat <- as.data.frame(fread(args$inspike))
  rownames(dat.spikeins.mat) <- dat.spikeins.mat$samp
  print("Peeking in dat spikeins mat")
  print(head(dat.spikeins.mat))
} else {
  stop("Not RData, rds of txt")
}

# Load data ---------------------------------------------------------------

# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562"
# mats.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   inname <- paste0("K562_count_tables_50000.", jmark, ".G1_G2_S.rds")
#   inf <- file.path(indir, inname)
#   readRDS(inf)
# })

mat <- readRDS(args$infile)


Y <- as.matrix(mat)
# Y <- as.matrix(mats.lst$H3K4me3)
# remove rows that are all zeros
rows.keep = which(rowSums(Y) > 0)

print(paste("Removing empty rows:", nrow(Y) - length(rows.keep)))

Y <- Y[which(rowSums(Y) > 0), ]

spikeincounts.sub <- dat.spikeins.mat[colnames(Y), ]

spikeincounts <- spikeincounts.sub$spikeincounts; names(spikeincounts) <- spikeincounts.sub$samp

assertthat::assert_that(ncol(Y) == length(spikeincounts))

print("Range:")
print(range(spikeincounts))
print(head(spikeincounts))

if (args$topndevgenes > nrow(Y)){
  print(paste("topnevgenes > nrow(Y), skipping deviance filter"))
  args$topndevgenes <- 0
}

print("Dim of Y before:")
print(dim(Y))
if (args$topndevgenes > 0){
  print(paste("Filtering topn genes by deviance:", args$topndevgenes))
  
  nvec <- spikeincounts.sub$spikeincounts 
  gdevs <- apply(mat, 1, function(xvec){
    scchicFuncs::binomial_deviance(x = xvec, p = sum(xvec) / sum(nvec), n = nvec)
  })
  genes.keep.vec <- sort(gdevs, decreasing = TRUE)[1:args$topndevgenes]
  genes.keep <- names(genes.keep.vec)
  Y <- Y[genes.keep, ]
} else {
  print("Taking all genes for glmpca")
}
print("Dim of Y after:")
print(dim(Y))

if (args$by_plate){
    if (!"plate" %in% colnames(spikeincounts.sub)){
        spikeincounts.sub <- spikeincounts.sub %>%
            rowwise() %>%
            mutate(plate = ClipLast(samp, jsep = "_"))
    }
    print(unique(spikeincounts.sub$plate))
    if (length(unique(spikeincounts.sub$plate)) > 1){
        X <- model.matrix(object = ~ plate, data = spikeincounts.sub)
    }  else {
        print("Number of factors in 'plate' is only 1, setting to NULL")
        X <- NULL
    }
} else {
    X <- NULL
}

print("Dimension of X")
print(dim(X))

system.time(
  glmpcaout <- glmpca::glmpca(Y = Y, L = args$K, fam = "poi", sz = spikeincounts, X = X, minibatch = args$minibatch, optimizer = args$optimizer, 
                              ctl = list(penalty = args$penalty, maxIter = args$maxIter, tol = args$tol))
)

save(glmpcaout, Y, spikeincounts, file = args$outfile)
