# Jake Yeung
# Date of Creation: 2021-02-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/9-check_hox_clusters.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_fits_all_bins_annotated_for_sup"


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




# Load hox clusters -------------------------------------------------------

inf.hox <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_rstudio/pdfs/topic_modules_50kb/topic24_H3K7me3_allmerged.2021-02-23.txt"
dat.hox <- fread(inf.hox)
dat.hox.filt <- subset(dat.hox, rnk < 100)

jmark.ref <- "H3K27me3"
params.wide.filt <- params.wide.lst[[jmark.ref]]
rows.keep <- rownames(params.wide.filt) %in% dat.hox.filt$bin

dat.params.filt.long <- data.frame(bin = rownames(params.wide.filt[rows.keep, ]), params.wide.filt[rows.keep, ], stringsAsFactors = FALSE) %>%
  melt() %>%
  dplyr::rename(logfc = value,
                param = variable) %>%
  rowwise() %>%
  mutate(ctype = gsub("Cluster", "", param),
         ctype = gsub(".Estimate", "", ctype))

ggplot(dat.params.filt.long, aes(x = ctype, y = logfc)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

jbin <- "chr15:7200000-7250000"
subset(dat.hox.filt, bin == jbin)

# Check single cells?  ----------------------------------------------------





