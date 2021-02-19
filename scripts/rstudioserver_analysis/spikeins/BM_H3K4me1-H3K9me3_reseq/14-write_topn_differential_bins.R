# Jake Yeung
# Date of Creation: 2021-02-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/14-write_topn_differential_bins.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(forcats)

library(topicmodels)
library(JFuncs)
library(scchicFuncs)

options(scipen=0)

make.plots <- TRUE

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"


# Constants ---------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load MATs  --------------------------------------------------------------

ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
ctypes.k9me3 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")

dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt")
  fread(inf)
})

dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})

jmetas.pretty.lst <- lapply(jmarks, function(jmark){
  jmeta <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes)
  } else { 
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes.k9me3)
  }
  jmeta <- jmeta %>% arrange(cluster, jrep)
})

cells.keep.lst <- lapply(jmetas.pretty.lst, function(jdat){
  jdat$cell
})


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
    bind_rows() %>%
    arrange(pval)
})


keepn <- 6085

bins.filt.lst <- lapply(pvals.lst, function(pvals.dat){
  pvals.dat[1:keepn, ]$bin
})

bins.bed.lst <- lapply(bins.filt.lst, function(jbins){
  GetBedFromCoords(jbins, add.chr = FALSE, strip.chr = TRUE)
})

# write to output
for (jmark in jmarks){
  print(jmark)
  jtmp <- bins.bed.lst[[jmark]]
  print(jtmp)
  fnametmp <- paste0("DE_bins_all_marks_top_", keepn, "_sorted_by_pval.", jmark, ".", Sys.Date(), ".bed")
  outftmp <- file.path(outdir, fnametmp)
  fwrite(jtmp, file = outftmp, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

