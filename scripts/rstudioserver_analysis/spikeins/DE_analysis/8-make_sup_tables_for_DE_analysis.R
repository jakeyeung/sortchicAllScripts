# Jake Yeung
# Date of Creation: 2021-02-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/8-make_sup_tables_for_DE_analysis.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


hubprefix <- "/home/jyeung/hub_oudenaarden"



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_fits_all_bins_annotated_for_sup"


# Get DE outputs ----------------------------------------------------------

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})


params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>% 
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  # make params more readable
  params.dat.all$ctype <- params.dat.all$param
  params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
  params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)
  return(params.dat.all)
})


pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, padj = p.adjust(xvec), stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})


# Make wide ---------------------------------------------------------------


params.wide.lst <- lapply(jmarks, function(jmark){
  params.wide <- as.data.frame(data.table::dcast(params.lst[[jmark]], formula = bin ~ param, value.var = "estimate"))
  params.wide$ClusterHSPCs.Estimate <- 0
  rownames(params.wide) <- params.wide$bin
  params.wide$bin <- NULL
  return(params.wide)
})

# add pvalhashes
pval.hash.lst <- lapply(jmarks, function(jmark){
  jhash <- hash::hash(pvals.lst[[jmark]]$bin, pvals.lst[[jmark]]$pval)
})



# Write outputs -----------------------------------------------------------

# filter DE bins 

infs.de.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

infs.high.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.de.bins.lst <- lapply(infs.de.lst, function(jinf){
  fread(jinf)
})

dat.high.bins.lst <- lapply(infs.high.lst, function(jinf){
  fread(jinf)
})

de.bins.lst <- lapply(dat.de.bins.lst, function(jdat){
  jdat$CoordOriginal
})

high.bins.lst <- lapply(dat.high.bins.lst, function(jdat){
  jdat$CoordOriginal
})


# Annotate bins  ----------------------------------------------------------


params.wide.annot.lst <- lapply(jmarks, function(jmark){
  jparam <- params.wide.lst[[jmark]]
  bnames <- rownames(jparam)
  pvals <- sapply(bnames, function(x) pval.hash.lst[[jmark]][[x]])
  jparam$pval <- pvals
  is.de <- bnames %in% de.bins.lst[[jmark]]
  is.high <- bnames %in% high.bins.lst[[jmark]]
  jparam$is.de <- is.de
  jparam$is.high <- is.high
  return(jparam)
})

# check
assertthat::assert_that(identical(lapply(de.bins.lst, length), lapply(params.wide.annot.lst, function(x) nrow(subset(x, is.de)))))
assertthat::assert_that(identical(lapply(high.bins.lst, length), lapply(params.wide.annot.lst, function(x) nrow(subset(x, is.high)))))

# write outputs
for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("SupplementalTable_PoissonFits_Outputs_DE_High_Annotated.", jmark, ".", Sys.Date(), ".txt"))
  write.table(params.wide.annot.lst[[jmark]], file = outf, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
}


# 
# 
# pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
# k9.bins <- which(pvals.lst2 < pvalcutoff)
# 
# pvals.adj.lst <- lapply(pvals.lst, function(jdat){
#   jdat$padj <- p.adjust(jdat$pval, method = "BH")
#   return(jdat)
# })
# 
# pfilt <- lapply(pvals.adj.lst, function(jdat){
#   subset(jdat, bin %in% names(k9.bins))
# })
# 

