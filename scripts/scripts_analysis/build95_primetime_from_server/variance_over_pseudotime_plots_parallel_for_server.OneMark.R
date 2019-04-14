# Jake Yeung
# Date of Creation: 2019-04-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_plots.R
# Plot the outputs 

rm(list=ls())

tstart <- Sys.time()

library(dplyr)
library(ggplot2)
library(tidyr)
library(JFuncs)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/FourierFunctions.R")
source("scripts/Rfunctions/VariabilityFunctions.R")

args <- commandArgs(trailingOnly=TRUE)
jmark <- args[[1]]
pdfout  <- args[[2]]

# Load constants ----------------------------------------------------------

jcol <- c("gray80", "gray50", "darkblue")
grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

jalpha <- 0.5

pseudo <- 0
jscale <- 1

mdpt.sd <- 1
lims.sd <- c(0, 3)
mdpt.fc <- 0.75
lims.fc <- c(0, 3)

jsize.facet <- 0.2
gw.jsize.facet <- 2

# do.plots <- TRUE

jstep <- 20000
# jtype <- "correlation"

pos.max <- 50 * 10^6
jstep <- 20000
lagmax <- pos.max / jstep

ci.interval <- 1.96  # corresponds to zscore 1.96 = alpha = 0.05 = 95% confidence interval

# Load data ---------------------------------------------------------------

# inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.RData"
# inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"
# inf <- "TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"

load(inf, v=T)

assertthat::assert_that(nrow(dat.umap.long.trajs[[jmark]]) > 0)
trajname <- ifelse(jmark == "H3K4me1", "bcell", "myeloid")
assertthat::assert_that(nrow(trajs[[jmark]][[trajname]]) > 0)

out <- MakeVariabilityPlots(jmark, trajname, tm.result.lst, dat.umap.long.trajs, pdfout = pdfout)

print(Sys.time() - tstart)
