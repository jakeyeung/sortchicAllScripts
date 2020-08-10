# Jake Yeung
# Date of Creation: 2020-06-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_poisson_model_gene_by_gene.R
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


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jdate <- "2020-06-09"

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir))

# dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.likeBM"))

outdir <- indir
outfits <- file.path(outdir, paste0("fit_poisson_model_on_TSS_ZF.", Sys.Date(), ".RData"))

infrdata <- paste0(jprefix, ".RData")

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)

dat.annots.filt.forfit <- lapply(dat.annot.lst.WKM, function(jdat){
  jdat <- subset(jdat, select = c(cell, cluster.new)) %>%
    mutate(Cluster = ifelse(cluster.new == "HSPCs", "aHSPCs", as.character(cluster.new)))  # set HSPC as intercept
  return(jdat)
})



print("Fitting... ")

system.time(
  jfits.lst.bymark <- lapply(jmarks, function(jmark){
    print(jmark)
    jmat.mark <- tss.mats.sc.filt.zf[[jmark]]
    dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
    ncuts.cells.mark <- ncuts.cells[[jmark]]
    cnames <- colnames(jmat.mark)
    jfits.lst <- apply(jmat.mark, MARGIN = 1, FUN = function(jrow){
      jout <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj = FALSE)
      return(jout)
    })
    return(jfits.lst)
  })
)

save(tss.mats.sc.filt.zf, dat.annots.filt.forfit, ncuts.cells, jfits.lst.bymark, file = outfits)

print(Sys.time() - jstart)



