# Jake Yeung
# Date of Creation: 2022-02-05
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/10-fit_pseudotime_glm_parallel_k4me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmark <- "k4me3"

FitGlmTrajs <- function(jrow, cnames, dat.annots, totalcounts.cells, returnobj = FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?
  
  # pvalue by deviance goodness of fit: https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  # offset is in log because the model says the log counts is equal to RHS
  # replace celltype-specific fit for m1
  
  if (!is.null(nrow(jrow))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(jrow), "rows"))
    row <- Matrix::colSums(jrow)
  }
  dat <- data.frame(cell = cnames, counts = jrow, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots, by = "cell") %>%
    left_join(., totalcounts.cells, by = "cell")
  
  m0.pois <- glm(counts ~ 1 + offset(log(totalcounts)), data = dat, family = "poisson")
  m1.pois <- glm(counts ~ 1 + ptime + offset(log(totalcounts)), data = dat, family = "poisson")
  
  # estimates2 <- summary(glm.complex)$coefficients[, "Estimate"]
  # names(estimates2) <- make.names(paste(names(estimates2), ".Estimate", sep = ""))
  # stderrors2 <- summary(glm.complex)$coefficients[, "Std. Error"]
  # names(stderrors2) <- make.names(paste(names(stderrors2), ".StdError", sep = ""))
  # 
  # dat.coefs.complex <- data.frame(estimate = estimates2[2:length(estimates2)],
  #                                 se = stderrors2[2:length(stderrors2)],
  #                                 intercept = estimates2[[1]],
  #                                 intercept.se = stderrors2[[1]])
  # dat.coefs.complex$param <- rownames(dat.coefs.complex)
  
  
  if (!returnobj){
    janova <- anova(m0.pois, m1.pois, test = "Chisq")
    # out.dat <- data.frame(pval10 = janova$`Pr(>Chi)`[[2]],
    #                       dev.diff10 = janova$`Deviance`[[2]],
    #                       df.diff10 = janova$`Df`[[2]],
    #                       t(as.data.frame(coefficients(m1.pois))),
    #                       stringsAsFactors = FALSE)
    estimates1 <- summary(m1.pois)$coefficients[, "Estimate"]
    names(estimates1) <- make.names(paste(names(estimates1), ".Estimate", sep = ""))
    stderrors1 <- summary(m1.pois)$coefficients[, "Std. Error"]
    names(stderrors1) <- make.names(paste(names(stderrors1), ".StdError", sep = ""))
    
    bic0 <- BIC(m0.pois)
    bic1 <- BIC(m1.pois)
    
    dat.coefs1 <- data.frame(estimate = estimates1[2:length(estimates1)],
                             se = stderrors1[2:length(stderrors1)],
                             intercept = estimates1[[1]],
                             intercept.se = stderrors1[[1]], 
                             pval10 = janova$`Pr(>Chi)`[[2]], 
                             dev.diff10 = janova$`Deviance`[[2]],
                             df.diff10 = janova$`Df`[[2]],
                             stringsAsFactors = FALSE)
    dat.coefs1$param <- rownames(dat.coefs1)
    
    out.dat <- list(bic0 = bic0,
                    bic1 = bic1,
                    coefs1 = dat.coefs1)
    return(out.dat)
  } else {
    return(list(m0.pois = m0.pois, m1.pois = m1.pois))
  }
}




# Load raw counts ---------------------------------------------------------


inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
load(inf.glmpca, v=T)
count.mat <- glm.inits$Y.filt


# Load metadata -----------------------------------------------------------

inf.traj.lst <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/trajs_outputs.", jmark, ".2022-02-05.rds")
traj.lst <- readRDS(inf.traj.lst)


# Fit pseudotime  ---------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits"

jctypes <- names(traj.lst); names(jctypes) <- jctypes
ncores <- 4

# jctype <- "Granulocytes"
# for (jctype in jctypes){
fit.out.lst.lst <- parallel::mclapply(jctypes, function(jctype){
  print(jctype)
  outf <- file.path(outdir, paste0("ptime_glm_fits_parallel.", jctype, ".", jmark, ".", Sys.Date(), ".rds"))
  dat.annots <- traj.lst[[jctype]] %>%
    filter(!is.na(ptime)) %>%
    dplyr::select(cell, ctype.from.LL, ptime)
  cells.keep <- dat.annots$cell
  count.mat.filt <- count.mat[, cells.keep]
  totalcounts.cells <- data.frame(cell = colnames(count.mat.filt), totalcounts = colSums(count.mat.filt), stringsAsFactors = FALSE)
  
  jgenes.vec <- rownames(count.mat.filt); names(jgenes.vec) <- jgenes.vec
  system.time(
    # fits.out.lst <- parallel::mclapply(jgenes.vec, function(jgene){
    fits.out.lst <- lapply(jgenes.vec, function(jgene){
      jrow <- count.mat.filt[jgene, ]
      cnames <- colnames(count.mat.filt)
      fit.out <- FitGlmTrajs(jrow = jrow, cnames = cnames, dat.annots = dat.annots, totalcounts.cells = totalcounts.cells, returnobj = FALSE)
      return(fit.out)
    # }, mc.cores = ncores)
    })
  )
  saveRDS(fits.out.lst, file = outf)
  print("Done fitting genes.")
  return(fits.out.lst)
}, mc.cores = length(jctypes))

outf2 <- file.path(outdir, paste0("final_ptime_glm_fits_parallel_allctypes.", jmark, ".", Sys.Date(), ".rds"))
saveRDS(fit.out.lst.lst, file = outf2)


