# Jake Yeung
# Date of Creation: 2020-06-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/downstream_poisson_model_analysis.SimulateCI.WrangleDownstream.R
# Wrangle CI into a manageable object to be used downstream

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# rename clusters
cnames.orig <- c("Clustereryth", "Clustergranu", "Clusterlymph")
cnames.new <- c("ClusterErythroblasts", "ClusterGranulocytes", "ClusterBcells")
cnames.hash <- hash::hash(cnames.orig, cnames.new)

# Load data  --------------------------------------------------------------

inf.ci <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits/fit_poisson_model_on_TSS_ZF.DownstreamWrangled.2020-06-10.ClusterRenamed.ConfidenceIntervals.smaller.RData"
outf.ci <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits/fit_poisson_model_on_TSS.ZF.DownstreamWrangled.2020-06-10.ConfidenceIntervals.smaller.Wrangled.RData"
assertthat::assert_that(!file.exists(outf.ci))

load(inf.ci, v=T)


# Wrangle  ----------------------------------------------------------------

jfits.fc.ci.long <- lapply(jmarks, function(jmark){
  rnames <- names(jfits.ci.lst[[jmark]])
  names(rnames) <- rnames
  jfits.tmp <- lapply(rnames, function(rname){
    jtmp <- jfits.ci.lst[[jmark]][[rname]]$params.fc.dat
    jtmp$bin <- rname
    return(jtmp)
  })  %>%
    bind_rows()
  jfits.tmp$mark <- jmark
  return(jfits.tmp)
}) %>%
  bind_rows()


jfits.int.ci <- lapply(jmarks, function(jmark){
  rnames <- names(jfits.ci.lst[[jmark]])
  names(rnames) <- rnames
  jfits.tmp <- lapply(rnames, function(rname){
    jtmp <- jfits.ci.lst[[jmark]][[rname]]$params.int.dat
    jtmp$bin <- rname
    return(jtmp)
  })  %>%
    bind_rows()
  jfits.tmp$mark <- jmark
  return(jfits.tmp)
}) %>%
  bind_rows()

# make mat of jfits.ci

jclsts <- unique(as.character(jfits.fc.ci.long$param))
names(jclsts) <- jclsts

jmat.fc.ci.lst <- lapply(jclsts, function(jclst){
  jmat.fc.lower <- reshape2::dcast(data = jfits.fc.ci.long %>% filter(param == jclst), formula = "bin ~ mark", value.var = "logLambdaLower")
  jmat.fc.upper <- reshape2::dcast(data = jfits.fc.ci.long %>% filter(param == jclst), formula = "bin ~ mark", value.var = "logLambdaUpper")
  colnames(jmat.fc.lower) <- c("bin", "H3K27me3.fc.lower", "H3K4me1.fc.lower", "H3K4me3.fc.lower")
  colnames(jmat.fc.upper) <- c("bin", "H3K27me3.fc.upper", "H3K4me1.fc.upper", "H3K4me3.fc.upper")
  jmat.fc.ci <- left_join(jmat.fc.lower, jmat.fc.upper, by = "bin")
})

names(jmat.fc.ci.lst) <- sapply(names(jmat.fc.ci.lst), function(x) AssignHash(x, cnames.hash, null.fill = x))


jmat.int.ci.lower <- reshape2::dcast(data = jfits.int.ci, formula = "bin ~ mark", value.var = "logLambdaLower")
jmat.int.ci.upper <- reshape2::dcast(data = jfits.int.ci, formula = "bin ~ mark", value.var = "logLambdaUpper")
colnames(jmat.int.ci.lower) <- c("bin", "H3K27me3.int.lower", "H3K4me1.int.lower", "H3K4me3.int.lower")
colnames(jmat.int.ci.upper) <- c("bin", "H3K27me3.int.upper", "H3K4me1.int.upper", "H3K4me3.int.upper")
jmat.int.ci <- left_join(jmat.int.ci.lower, jmat.int.ci.upper, by = "bin")

# Save outputs ------------------------------------------------------------

save(jmat.fc.ci.lst, jmat.int.ci, file = outf.ci)
