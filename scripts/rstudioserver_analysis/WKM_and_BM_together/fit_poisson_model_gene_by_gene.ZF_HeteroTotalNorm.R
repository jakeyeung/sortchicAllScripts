# Jake Yeung
# Date of Creation: 2020-06-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.ZF_HeteroTotalNorm.R
# description


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

hubprefix <- "/home/jyeung/hub_oudenaarden"

jdists <- c(50000L, 10000L)
jnorms <- c("ByTotalFromBins", "ByHetero")


# jdist <- "10000"
# # jdist <- "50000"
# jnorm <- "ByHetero"
# # jnorm <- "ByTotalFromBins"

jcell <- "PZ-ChIC-ZFWKM-H3K4me1-G1-4_373"

for (jdist in jdists){
  print(jdist)
  jfactor <- jdist / 10000L  # factor to multiply to ncuts total in order to make different winsizes comparable
  for (jnorm in jnorms){
    
    print(jnorm)
    jspecies <- "ZebrafishWKM"
    jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
    outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists"
    dir.create(outdir)
    outfits <- file.path(outdir, paste0("fit_poisson_model_on_TSS.", jspecies, ".NormMeth_", jnorm, ".dist_", jdist, ".RData"))
    
    # assertthat::assert_that(!file.exists(outfits))
    if (file.exists(outfits)){
      print(paste("Skipping:", jdist, jnorm))
      next
    }
    
    inf.offsets <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts/", jspecies, "_HeteroTotalCounts_50kb_bins.2020-06-16.RData")
    
    
    # Load Annots and TSS transcripts I keep ----------------------------------------------------
    
    jdate <- "2020-06-09"
    indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
    assertthat::assert_that(dir.exists(indir))
    # dir.create(indir, showWarnings = TRUE)
    jprefix <- file.path(indir, paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.likeBM"))
    infrdata <- paste0(jprefix, ".RData")
    load(infrdata, v=T)
    
    # take transcripts: can be different window siezs
    refmark <- "H3K4me3"
    rnames <- rownames(tss.mats.sc.filt.zf[[refmark]])
    tx.keep <- sapply(rnames, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)
    
    
    dat.annots.filt.forfit <- lapply(dat.annot.lst.WKM, function(jdat){
      # jdat$cluster.new <- as.character(jdat$cluster.new)
      jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
        mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", as.character(cluster.new)))  # set HSPC as intercept
      return(jdat)
    })
    
    cells.clstrfilt <- lapply(dat.annots.filt.forfit, function(jdat){
      return(jdat$cell)
    })
    
    load(inf.offsets, v=T)
    
    if (jnorm == "ByHetero"){
      ncuts.for.fit <- lapply(dat.ncuts.hetero.total, function(jdat){
        # expects ncuts.total
        subset(jdat, select = c(cell, ncuts.hetero)) %>%
          dplyr::rename(ncuts.total = ncuts.hetero)
      })
    } else if (jnorm == "ByTotalFromBins"){
      ncuts.for.fit <- lapply(dat.ncuts.hetero.total, function(jdat){
        # expects ncuts.total
        subset(jdat, select = c(cell, ncuts.total))
      })
    } else {
      warning("jnorm must be ByHetero or ByTotalFromBins", jnorm)
    }
    
    # adjust by jfactor
    ncuts.for.fit <- lapply(ncuts.for.fit, function(jdat){
      jdat$ncuts.total * jfactor  # multiply by factor so different winsizes are comparable
      return(jdat)
    })
 
    
    
    # Load the new tables ------------------------------------------------
    
    # must be TSS tables
    indir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.CodingOnly.imputevarfilt.lessstringent.mapq_40.winsize_", jdist, ".r1only")
    tss.mats.filt.fromref.cellfilt <- lapply(jmarks, function(jmark){
      print(jmark)
      # fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.csv")
      fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv")
      inf <- file.path(indir, fname)
      mat <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE)
      # filter rownames to match tx.keep
      tx.all <- sapply(rownames(mat), function(rname) paste(strsplit(rname, ";")[[1]][2:3], collapse = ";"))
      tx.filt <- tx.all %in% tx.keep
      return(mat[tx.filt, cells.clstrfilt[[jmark]]])
    })
    
    # rcommon <- Reduce(intersect, lapply(jmats, rownames))
    # tss.mats.filt.fromref.cellfilt <- lapply(jmats, function(jmat) jmat[rcommon, ])
    
    lapply(tss.mats.filt.fromref.cellfilt, dim)
    
    
    # Fit data ----------------------------------------------------------------
    
    lapply(ncuts.for.fit, function(jdat) range(jdat$ncuts.total))
    
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
          # jfits.lst <- lapply(jrow.names, function(jrow.name){
          # print(jrow.name)
          jrow <- jmat.mark[jrow.name, ]
          jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = FALSE)
          # jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.for.fit.mark, jbin = NULL, returnobj = TRUE)  # for debugging
          return(jout)
        }, mc.cores = ncores)
        # })
        return(jfits.lst)
      })
    )
    
    print("Saving objects")
    save(tss.mats.filt.fromref.cellfilt, dat.annots.filt.forfit, ncuts.for.fit, jfits.lst.bymark, file = outfits)
    print("Done Saving objects")
    
  }
}



print(Sys.time() - jstart)



