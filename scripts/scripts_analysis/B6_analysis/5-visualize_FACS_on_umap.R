# Jake Yeung
# Date of Creation: 2019-05-10
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/5-visualize_FACS_on_umap.R
# Visualize FACS on umap 

rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(hash)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxB6.R")

# Functions ---------------------------------------------------------------




# Constants ---------------------------------------------------------------

# jmark <- "H3K9me3"
# jmark <- "H3K27me3"
jmark <- "H3K4me1"
# jmark <- "H3K4me3"

# Load umaps --------------------------------------------------------------

inf.dat <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.", jmark, ".RData")
assertthat::assert_that(file.exists(inf.dat))

load(inf.dat, v=T)

dat.umap.long$repl <- as.character(sapply(as.character(dat.umap.long$cell), function(x) strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
dat.umap.long$techname <- as.character(sapply(as.character(dat.umap.long$cell), function(x) strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
dat.umap.long$cellindx <- as.character(sapply(as.character(dat.umap.long$cell), function(x) strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[2]]))
rep2techname <- hash(dat.umap.long$repl, dat.umap.long$techname)


# Load reference facs to match colnames -----------------------------------

inf.facs.ref <- "/Users/yeung/data/scchic/facs/H3K9me3_index_2mice_4plates.csv"
dat.refs <- read.table(inf.facs.ref)
cnames.noX <- colnames(dat.refs)[!grepl(pattern = "^X", colnames(dat.refs))]


# Load FACS ---------------------------------------------------------------

prefix <- paste0("B6-13W1-BM-", jmark, "-")

if (jmark != "H3K4me1"){
  repl.vec <- seq(4)
} else {
  repl.vec <- c(2, 3, 4)
}

dat.facs.filt <- lapply(repl.vec, function(repl) LoadFACS(jmark, repl, cnames.noX, rep2techname)) %>%
  bind_rows()
rownames(dat.facs.filt) <- dat.facs.filt$cellname

dat.pca <- GetFACSloadings(dat.facs.filt %>% dplyr::select(-cell), PC = seq(5))

dat.facs.filt$loadings <- GetFACSloadings(dat.facs.filt %>% dplyr::select(-cell), PC = 1, make.pos = TRUE)

dat.merge <- left_join(dat.umap.long, dat.facs.filt %>% dplyr::select(cell, loadings))

PlotXYWithColor(dat.merge %>% filter(!is.na(loadings)), xvar = "umap1", yvar = "umap2", cname = "loadings")



# Notes -------------------------------------------------------------------


# repl <- 1
# inf.facs <- paste0("/Users/yeung/data/scchic/facs/20190403/20190403_Mouse_BM_B6_", jmark, "_00", repl, "_index.csv")
# dat.facs <- read.csv(inf.facs)
# dat.facs$X <- NULL  # confuses downstream colnames
# dat.facs$repl <- as.character(repl)
# dat.facs$techname <- sapply(as.character(dat.facs$repl), function(x) rep2techname[[x]])
# dat.facs$mark <- jmark
# dat.facs$wellindx <- as.numeric(c(seq(1:356), seq(360:379)+360))
# dat.facs$wellindx0 <- dat.facs$wellindx - 1  # zero based
# dat.facs$cell <- paste0(prefix, dat.facs$repl, "_", dat.facs$techname, "_", dat.facs$wellindx0)
# # keep colnames that start with X or in cnamesnoX
# cnames.keep <- colnames(dat.facs)[grepl(pattern = "^X", colnames(dat.facs))]
# # add noX
# cnames.keep <- c(cnames.noX, cnames.keep)
# 
# dat.facs.filt <- dat.facs %>% dplyr::select(cnames.keep)
# rownames(dat.facs.filt) <- dat.facs$cell
# 
# X <- prcomp(dat.facs.filt, center = TRUE, scale. = TRUE)
# 
# plot(X$x[, 1], X$x[, 2])
# 
# dat.facs$loading <- X$x[, 1]

