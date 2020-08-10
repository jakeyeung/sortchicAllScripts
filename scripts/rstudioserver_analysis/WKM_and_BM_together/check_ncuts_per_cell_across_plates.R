# Jake Yeung
# Date of Creation: 2020-06-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/check_ncuts_per_cell_across_plates.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
fewer.k27me3 <- TRUE

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", "2020-06-05"))

outdir <- indir
outfits <- file.path(outdir, paste0("fit_poisson_model_on_TSS.RData"))

infrdata <- paste0(jprefix, ".smaller.RData")

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)

ncuts.cells.plates <- lapply(jmarks, function(jmark){
  ncuts <- ncuts.cells[[jmark]]
  ncuts$plate <- sapply(ncuts$cell, ClipLast, jsep = "_")
  ncuts$cond <- sapply(ncuts$cell, GetCondFromSamp, mark = jmark)
  return(ncuts)
})


# Plot across plates ------------------------------------------------------

m.lst <- lapply(jmarks, function(jmark){
  jdat <- ncuts.cells.plates[[jmark]]
  m <- ggplot(jdat, aes(x = forcats::fct_reorder(.f = plate, .x = log10(ncuts.total), .fun = median, .desc = TRUE), y = log10(ncuts.total), fill = cond)) + 
    geom_boxplot() + 
    # geom_violin() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +   
    xlab("") + 
    ggtitle(jmark)
  return(m)
})
print(m.lst)

multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)
