# Jake Yeung
# Date of Creation: 2020-11-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/10-de_genes_for_Buys.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(scchicFuncs)

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

e2g.hash <- hash::invert(g2e.hash)

dat.sorted.norm.long$genesymbol <- sapply(dat.sorted.norm.long$gene, function(x) AssignHash(x = x, jhash = e2g.hash, null.fill = NA))

dat.sorted.norm.long.filt <- subset(dat.sorted.norm.long, !is.na(genesymbol))

# make table for Buys

jannots <- names(de.ens.sorted.stringent)
names(jannots) <- jannots

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/gene_lists"
# lapply(jannots, function(jannot){
for (jannot in jannots){
  print(jannot)
  fname <- paste0("DE_genes_BM.", jannot, ".txt")
  jgenes <- as.character(de.ens.sorted.stringent[[jannot]])
  jsub <- subset(dat.sorted.norm.long.filt, gene %in% jgenes)
  print(dim(jsub))
  fwrite(jsub, file = file.path(outdir, fname), sep = "\t", quote = FALSE)
}



# 
# infrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression/integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_TRUE.forPoissonRegression.CountR1only.2020-06-05.smaller.RData"
# assertthat::assert_that(file.exists(infrdata))
# load(infrdata, v=T)
# 
# jannots <- names(de.ens.sorted.stringent)
# names(jannots) <- jannots