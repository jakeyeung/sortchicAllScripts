# Jake Yeung
# Date of Creation: 2021-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/1-check_top_6085_bins.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Check 6085 bins first time ----------------------------------------------

bins1.lst <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.", jmark, ".2021-02-15.bed")
  fread(inf)
})

bins2.lst <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_top_6085_sorted_by_pval.", jmark, ".2021-02-19.bed")
  fread(inf)
})

common.bins.lst <- lapply(jmarks, function(jmark){
  common.bins <- intersect(bins1.lst$H3K9me3$V4, bins2.lst$H3K9me3$V4)
})


# Check pvalues from fits -------------------------------------------------


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

pvals.adj.lst <- lapply(pvals.lst, function(jdat){
  jdat$padj <- p.adjust(jdat$pval, method = "BH")
  return(jdat)
})



# Check pvalues for top bins  ---------------------------------------------

pvals.range.lst <- lapply(jmarks, function(jmark){
  jbins <- bins1.lst[[jmark]]$V4
  pvals.range <- range(subset(pvals.adj.lst[[jmark]], bin %in% jbins)$pval)
})
