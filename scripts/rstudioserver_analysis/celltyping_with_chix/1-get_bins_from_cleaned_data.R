# Jake Yeung
# Date of Creation: 2020-08-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/celltyping_with_chix/1-get_bins_from_cleaned_data.R
# Get bins from cleaned data 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load data  --------------------------------------------------------------

jmark <- "H3K4me1"
hubpath <- "/home/jyeung/hub_oudenaarden"
outdir <- file.path(hubpath, "jyeung/data/scChiC/bed_from_count_mat")
outf <- file.path(outdir, paste0(jmark, "_peaks_from_hiddendomains.bed"))
outf.nochromo <- file.path(outdir, paste0(jmark, "_peaks_from_hiddendomains.nochromo.bed"))


inf.lda <- file.path(hubpath, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains.NewCountFilters/lda_outputs.merged.H3K27me3.minlength_1000.cutoff_analysis.merged.withchr.annotated.NewCountFilters.K-30.binarize.FALSE/ldaOut.merged.H3K27me3.minlength_1000.cutoff_analysis.merged.withchr.annotated.NewCountFilters.K-30.Robj")

load(inf.lda, v=T)

bins.keep <- rownames(count.mat)
coords.keep <- sapply(bins.keep, function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE)


# Write bedfile  ----------------------------------------------------------

dat.bins <- data.frame(chromo = sapply(coords.keep, GetChromo),
                       jstart = as.numeric(sapply(coords.keep, GetStart)), 
                       jend = as.numeric(sapply(coords.keep, GetEnd)), 
                       jname = bins.keep, 
                       stringsAsFactors = FALSE)

dat.bins.nochromo <- dat.bins %>%
  ungroup() %>%
  mutate(chromo = gsub("^chr", "", chromo))

fwrite(dat.bins, file = outf, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(dat.bins.nochromo, file = outf.nochromo, sep = "\t", row.names = FALSE, col.names = FALSE)

