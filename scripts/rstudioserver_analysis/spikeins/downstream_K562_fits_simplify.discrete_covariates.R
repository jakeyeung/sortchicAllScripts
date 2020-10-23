# Jake Yeung
# Date of Creation: 2020-08-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/downstream_K562_fits.discrete_covariates.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(parallel)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"
jsuffix <- "topn_5000.glmpcaout.penalty_5.by_plate.RData"
jprefix <- ".G1_G2_S."

# jmark <- "H3K27me3"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks


mclapply(jmarks, function(jmark){

    print(jmark)

    # Load inputs -------------------------------------------------------------

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

    inf.fits <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.discrete_covariate/fit_cellcycle_pseudotime.", jmark, ".2020-08-16.RData"))
    load(inf.fits, v=T)

    outf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.discrete_covariate/fit_cellcycle_pseudotime.", jmark, ".2020-08-16.simplified.RData"))


    # Get slope  --------------------------------------------------------------


    slopes.S <- lapply(jfits, function(jfit){
      return(coef(jfit)[["cellcycle.str1_S"]])
    })

    slopes.G2 <- lapply(jfits, function(jfit){
      return(coef(jfit)[["cellcycle.str2_G2/M"]])
    })

    pvals.S <- lapply(jfits, function(jfit){
      return(summary(jfit)$coefficients["cellcycle.str1_S", "Pr(>|z|)"])
    })

    pvals.G2 <- lapply(jfits, function(jfit){
      return(summary(jfit)$coefficients["cellcycle.str2_G2/M", "Pr(>|z|)"])
    })

    coefs.dat <- data.frame(coord = names(slopes.S), 
                            slope.S = unlist(slopes.S), pval.S = unlist(pvals.S), 
                            slope.G2 = unlist(slopes.G2), pval.G2 = unlist(pvals.G2), 
                            stringsAsFactors = FALSE)

    save(coefs.dat, file = outf)

}, mc.cores = length(jmarks))

# for (jmark in jmarks){

# }

