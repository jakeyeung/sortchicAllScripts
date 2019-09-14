# Jake Yeung
# Date of Creation: 2019-08-11
# File: ~/projects/scchic/scripts/scripts_analysis/revisions/create_bin_bed.R
# Make bed file for bins used in analysis so we can compare it with the bins in the ChIP-seq.


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gganimate)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Load data ---------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/pdfs/revision_figures"

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
infs.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"

infs$H3K4me3 <- infs.stringent

tm.result.lst <- lapply(infs, LoadGetTmResult)

# Make bin list -----------------------------------------------------------

jbins.lst <- lapply(tm.result.lst, function(lst){
  colnames(lst$terms)
})

# check all elements are equal 
# https://stackoverflow.com/questions/34313463/test-for-equality-between-all-members-of-list
assertthat::assert_that(Reduce(function(x, y) x && identical(y, jbins.lst[[1]]), init = TRUE, jbins.lst))

# take first one because it doesn't matter
jbins <- jbins.lst[[1]]

# Write to file -----------------------------------------------------------

out.dat <- data.frame(chromo = sapply(jbins, GetChromo), 
                      start = sapply(jbins, GetStart), 
                      end = sapply(jbins, GetEnd),
                      stringsAsFactors = FALSE)

fwrite(out.dat, file = paste0("/Users/yeung/data/scchic/revision_objs/BM_bins.", Sys.Date(), ".bed"), quote = FALSE, sep = "\t", col.names = FALSE)


