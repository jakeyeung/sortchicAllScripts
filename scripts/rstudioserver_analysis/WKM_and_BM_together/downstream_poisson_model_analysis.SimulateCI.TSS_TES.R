# Jake Yeung
# Date of Creation: 2020-08-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/downstream_poisson_model_analysis.SimulateCI.TSS_TES.R
# 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(ggrastr)
library(scchicFuncs)

# Functions ---------------------------------------------------------------



# Constants ---------------------------------------------------------------

ncores <- 16
jdate <- "2020-06-05"

# jnorm <- "ByHetero"
# jnorm <- "ByTotalFromBins"

fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
assertthat::assert_that(dir.exists(indir))

jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", jdate))

# datdir <- "/home/jyeung/data/from_rstudioserver/poisson_fits"
# datdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm"
datdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists.TSS_TES"
dir.create(datdir)
outfits <- file.path(datdir, paste0("fit_poisson_model_on_TSS.MouseBM.TSS_TES.RData"))
# outfits <- file.path(datdir, paste0("fit_poisson_model_on_TSS.BM.NormMeth_", jnorm, ".RData"))
outfits.wrangled <- file.path(datdir, paste0("fit_poisson_model_on_TSS.MouseBM.TSS_TES.CI.RData"))
# outfits.wrangled <- file.path(datdir, paste0("fit_poisson_model_on_TSS.BM.NormMeth_", jnorm, ".CI.RData"))

assertthat::assert_that(file.exists(outfits))
assertthat::assert_that(!file.exists(outfits.wrangled))

infrdata <- paste0(jprefix, ".smaller.RData")

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)
load(outfits, v=T)

# pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)


# Wrangle -----------------------------------------------------------------

jfits.dat.lst <- lapply(jmarks, function(jmark){
  jfit.dat <- do.call(rbind, jfits.lst.bymark[[jmark]])
  jfit.dat$bin <- rownames(jfit.dat)
  jfit.dat$gene <- sapply(as.character(jfit.dat$bin), function(b) strsplit(strsplit(b, split = "\\.\\.")[[1]][[2]], split = "_")[[1]][[1]])
  jfit.dat$ens <- sapply(jfit.dat$gene, function(g) AssignHash(g, g2e, null.fill = g))
  return(jfit.dat)
})

jfits.long.lst <- lapply(jmarks, function(jmark){
  jlong <- jfits.dat.lst[[jmark]] %>%
    dplyr::select(-c(dev.diff, df.diff)) %>%
    reshape2::melt(., id.vars = c("bin", "gene", "ens", "pval"), variable.name = "cluster", value.name = "logLambda")
  jlong$mark <- jmark
  # extract the X intercept as separate column 
  jlong.noint <- subset(jlong, cluster != "X.Intercept.")
  jlong.int <- subset(jlong, cluster == "X.Intercept.")  %>%
    dplyr::rename(logintercept = logLambda) %>%
    dplyr::select(bin, logintercept, mark)
  jlong.merge <- left_join(jlong.noint, jlong.int)
  return(jlong.merge)
})

print(paste("Running fits... ncores:", ncores))


# # debug
# jmark <- "H3K27me3"
# jrow.name <- "chr4:44705408-44715408;Pax5;5"

system.time(
  jfits.ci.lst <- lapply(jmarks, function(jmark){
    jrow.names <- rownames(tss.mats.filt.fromref.cellfilt[[jmark]])
    names(jrow.names) <- jrow.names
    cnames <- colnames(tss.mats.filt.fromref.cellfilt[[jmark]])
    dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
    ncuts.for.fit.mark <- ncuts.for.fit[[jmark]]
    # ncuts.for.fit.mark <- ncuts.for.fit[[jmark]] %>%
    #   mutate(ncuts.total = 1000000000)
    # print(head(ncuts.for.fit.mark))
      
    refits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
      refit <- RefitPoissonForPlot(jrow = tss.mats.filt.fromref.cellfilt[[jmark]][jrow.name, ], cnames = cnames, dat.annots.filt.mark = dat.annots.filt.mark, 
                                   ncuts.cells.mark = ncuts.for.fit.mark, return.means = FALSE)
      return(refit)
    }, mc.cores = ncores)
  })
)

# if (!file.exists(outfits.wrangled)){
  save(jfits.ci.lst, file = outfits.wrangled)
# }

