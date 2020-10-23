# Jake Yeung
# Date of Creation: 2020-08-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/downstream_K562_fits.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"
jsuffix <- "topn_5000.glmpcaout.penalty_5.by_plate.RData"
jprefix <- ".G1_G2_S."

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

for (jmark in jmarks){
    print(jmark)

    # Load inputs -------------------------------------------------------------

    inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".2020-08-14.RData"))
    outf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle/fit_cellcycle_pseudotime.", jmark, ".2020-08-14.simplified.RData"))
    assertthat::assert_that(!file.exists(outf.fits))

    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, jprefix, jsuffix))
    assertthat::assert_that(file.exists(inf.glmpca))


    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/K562_count_tables_50000.", jmark, jprefix, "rds"))
    inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")

    assertthat::assert_that(file.exists(inf))
    assertthat::assert_that(file.exists(inf.spike))

    mat <- readRDS(inf)
    load(inf.spike, v=T)
    load(inf.glmpca, v=T)


    # Load fits ---------------------------------------------------------------


    load(inf.fits, v=T)


    # Get slope  --------------------------------------------------------------

    slopes <- lapply(jfits, function(jfit){
      return(coef(jfit)[["pseudotime"]])
    })

    pvals <- lapply(jfits, function(jfit){
      return(summary(jfit)$coefficients["pseudotime", "Pr(>|z|)"])
    })

    coefs.dat <- data.frame(coord = names(slopes), slope = unlist(slopes), pval = unlist(pvals), stringsAsFactors = FALSE)

    save(coefs.dat, file = outf.fits)

}

