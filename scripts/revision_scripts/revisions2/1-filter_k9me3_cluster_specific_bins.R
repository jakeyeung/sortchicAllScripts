# Jake Yeung
# Date of Creation: 2022-07-21
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/1-filter_k9me3_cluster_specific_bins.R
# copied largely from /nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/10-DE_downstream_all_marks.K9me3_bins.R
# Create cluster-specific bins

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# keeptop <- 150
keeptop <- 500

# Load bins ---------------------------------------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/H3K9me3_H3K4me1_analysis"

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(indir, paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
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
  return(params.dat.all)
})


pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < 1e-10)



# Get k9 bins and plot  ---------------------------------------------------


pval.k9.sub <- subset(pvals.lst$H3K9me3, pval < 1e-10) %>%
  arrange(desc(pval))

k9.bins.names <- unique(pval.k9.sub$bin)
ctypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", ctypes.keep, ".Estimate", sep = "")

params.dat.wide <- data.table::dcast(subset(params.lst$H3K9me3, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate") %>%
  rowwise() %>%
  mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
         ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
         ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
         ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
         Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
         Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
         Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
         HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))

params.dat.wide.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(params.lst[[jmark]], bin %in% k9.bins.names & param %in% params.keep) %>%
    group_by(bin) %>% filter(max(abs(estimate)) < 5)
  jdat <- data.table::dcast(jsub,
                            formula = bin ~ param, value.var = "estimate") %>%
    rowwise() %>%
    mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
           ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
           ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
           ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
           ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
           ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
           Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
           Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
           Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
           HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))
  # keep only effect cnames
  cnames.keep.i <- grep("effect$", colnames(jdat))
  cnames.new <- paste(colnames(jdat)[cnames.keep.i], jmark, sep = "_")
  colnames(jdat)[cnames.keep.i] <- cnames.new
  cnames.keep.bin.i <- grep("bin", colnames(jdat))
  cnames.keep.merged.i <- c(cnames.keep.bin.i, cnames.keep.i)
  jdat.filt <- jdat[, cnames.keep.merged.i]
  return(jdat.filt)
})



# Get clusterspecific -----------------------------------------------------

jsort.hspcs <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(HSPCs.effect)
  arrange(desc(HSPCs.effect))
jbins.hspcs <- jsort.hspcs$bin[1:keeptop]

jsort.bcell <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(desc(Bcells.effect))
  arrange(Bcells.effect)
jbins.bcell <- jsort.bcell$bin[1:keeptop]

jsort.granu <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(desc(Granulocytes.effect))
  arrange(Granulocytes.effect)
jbins.granu <- jsort.granu$bin[1:keeptop]

jsort.eryth <- params.dat.wide %>%
  group_by(bin) %>%
  # arrange(descEryths.effect))
  arrange(Eryths.effect)
jbins.eryth <- jsort.eryth$bin[1:keeptop]



# Downstream --------------------------------------------------------------


bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

bins.keep.lst <- list("Eryths" = jbins.eryth,
                      "Bcells" = jbins.bcell,
                      "Granulocytes" = jbins.granu,
                      "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

length(bins.keep)



# Write bed table ---------------------------------------------------------

chromos <- sapply(bins.keep, function(b) gsub("^chr", "", strsplit(b, ":")[[1]][[1]]))
startends <- sapply(bins.keep, function(b) strsplit(b, ":")[[1]][[2]])
starts <- sapply(startends, function(se) strsplit(se, "-")[[1]][[1]])
ends <- sapply(startends, function(se) strsplit(se, "-")[[1]][[2]])
bins.keep.nochromo <- gsub("^chr", "", bins.keep)

dat.bins <- data.frame(Chr = chromos, Start = starts, End = ends, Bin = bins.keep.nochromo, stringsAsFactors = FALSE)

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/H3K9me3_H3K4me1_analysis"
fwrite(dat.bins, file = file.path(outdir, paste0("H3K9me3_cluster_specific_bins.keeptop_", keeptop, ".", Sys.Date(), ".txt")), sep = "\t", col.names = FALSE)



# load("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/H3K9me3_H3K4me1_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2020-12-20.H3K27me3rep2rep3reseq.RData", v=T)



