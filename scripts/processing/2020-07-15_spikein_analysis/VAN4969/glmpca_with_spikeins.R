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

load(args$inspike, verbose=TRUE)

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
  glmpcaout <- glmpca::glmpca(Y = Y, L = args$K, fam = "poi", penalty = args$penalty, sz = spikeincounts, X = X)
)

save(glmpcaout, Y, spikeincounts, file = args$outfile)
