# Jake Yeung
# Date of Creation: 2020-06-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.BM_HeteroTotalNorm.R
# 



rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(JFuncs)
library(scchicFuncs)

ncores <- 24


# jdists <- c(50000L, 10000L)
jdists <- c(10000L)
# jnorms <- c("ByTotalFromBins", "ByHetero")
jnorms <- c("ncuts.inbins", "ncuts.alltss")


# jdist <- 10000
# jnorm <- "ByHetero"
# # jnorm <- "ByTotalFromBins"

for (jdist in jdists){
  jfactor <- jdist / 10000
  
  for (jnorm in jnorms){
    print(paste(jnorm, jdist))
    
    fewer.k27me3 <- TRUE
    jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
    
    indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
    dir.create(indir, showWarnings = TRUE)
    jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", "2020-06-05"))
    
    outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists"
    dir.create(outdir)
    outfits <- file.path(outdir, paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".bsize_", jdist, ".RData"))
    
    assertthat::assert_that(!file.exists(outfits))
    
    infrdata <- paste0(jprefix, ".smaller.RData")
    
    # inf.offsets <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts/MouseBM_HeteroTotalCounts_50kb_bins.2020-06-16.RData"
    # inf.offsets <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts/MouseBM_HeteroTotalCounts_50kb_bins.2020-06-16.RData"
    inf.offsets <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bins_tss_totalcuts.countsByCluster/MouseBMFromTopics.2000.bins_tss_hetero_and_totalcuts.celltable.txt"
    
    assertthat::assert_that(file.exists(infrdata))
    
    # load(inf.offsets, v=T)
    dat.offsets <- fread(inf.offsets)
    assertthat::assert_that(jnorm %in% colnames(dat.offsets))
    dat.offsets.lst <- split(x = dat.offsets, f = dat.offsets$mark)
    
    # split into list
    load(infrdata, v=T)
    
    cells.keep <- lapply(tss.mats.filt.fromref.cellfilt, function(jmat){
      colnames(jmat)
    })
    
    # take transcripts: can be different window siezs
    refmark <- "H3K4me3"
    rnames <- rownames(tss.mats.filt.fromref.cellfilt[[refmark]])
    tx.keep <- sapply(rnames, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)
    
    
    rm(tss.mats.filt.fromref.cellfilt)  # load a new one later
    
    
    dat.annots.filt.forfit <- lapply(dat.annots.filt, function(jdat){
      jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
        mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", as.character(cluster.new)))  # set HSPC as intercept
      return(jdat)
    })
    
    cells.clstrfilt <- lapply(dat.annots.filt.forfit, function(jdat){
      return(jdat$cell)
    })
    
    ncuts.for.fit <- lapply(dat.offsets.lst, function(jdat){
      # expects ncuts.total in model so rename
      cnames.keep <- c("cell", jnorm)
      cnames.output <- c("cell", "ncuts.total")
      cnames.keep.i <- which(colnames(jdat) %in% cnames.keep)
      jdat.select <- jdat[, ..cnames.keep.i] 
      colnames(jdat.select) <- cnames.output
      return(jdat.select)
    })
    
    
    # adjust by jfactor
    ncuts.for.fit <- lapply(ncuts.for.fit, function(jdat){
      jdat$ncuts.total * jfactor  # multiply by factor so different winsizes are comparable
      return(jdat)
    })
    
    
    # Load the new tables ------------------------------------------------
    
    # must be TSS tables
    indir.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/debug/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.all_tx.noR2.Buys"
    tss.mats.filt.fromref.cellfilt <- lapply(jmarks, function(jmark){
      print(jmark)
      # fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv")
      fname <- paste0(jmark, ".countTableTSS.mapq_40.TSS_", jdist, ".blfiltered.csv")
      inf <- file.path(indir.tss, fname)
      mat <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE)
      # filter rownames to match tx.keep
      tx.all <- sapply(rownames(mat), function(rname) paste(strsplit(rname, ";")[[1]][2:3], collapse = ";"))
      tx.filt <- tx.all %in% tx.keep
      return(mat[tx.filt, cells.clstrfilt[[jmark]]])
    })
    
    
    # Fit ---------------------------------------------------------------------
    
    
    
    print("Fitting... ")
    
    system.time(
      jfits.lst.bymark <- lapply(jmarks, function(jmark){
        # jfits.lst.bymark <- parallel::mclapply(jmarks, function(jmark){
        print(jmark)
        jmat.mark <- tss.mats.filt.fromref.cellfilt[[jmark]]
        dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
        ncuts.for.fit.mark <- ncuts.for.fit[[jmark]]
        cnames <- colnames(jmat.mark)
        
        jrow.names <- rownames(jmat.mark)
        names(jrow.names) <- jrow.names
        jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
          jrow <- jmat.mark[jrow.name, ]
          jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE)
          return(jout)
        }, mc.cores = ncores)
        
        # jfits.lst <- apply(jmat.mark, MARGIN = 1, FUN = function(jrow){
        #   jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE)
        #   return(jout)
        # })
        return(jfits.lst)
        # }, mc.cores = ncores)
      })
    )
    
    print("Saving objects")
    save(tss.mats.filt.fromref.cellfilt, dat.annots.filt.forfit, ncuts.for.fit, jfits.lst.bymark, file = outfits)
    print("Done Saving objects")
    print(Sys.time() - jstart)
    
  }
}




print(Sys.time() - jstart)



