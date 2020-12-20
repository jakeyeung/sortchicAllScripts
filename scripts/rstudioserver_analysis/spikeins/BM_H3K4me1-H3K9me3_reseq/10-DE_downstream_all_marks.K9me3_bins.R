# Jake Yeung
# Date of Creation: 2020-12-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/10-DE_downstream_all_marks.K9me3_bins.R
# Do DE downstream 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

keeptop <- 150
low.in.k9 <- TRUE
# low.in.k9 <- FALSE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


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

params.dat.wide <- data.table::dcast(subset(params.lst$H3K9me3, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate") %>%
  rowwise() %>%
  mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
         ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
         ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
         ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
         bcell.effect = ClusterBcells.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
         eryth.effect = ClusterEryths.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
         granu.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
         mean.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))


if (low.in.k9){
  jsort.hspcs <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(mean.effect)
    arrange(desc(mean.effect))
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
  
  jsort.bcell <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(desc(bcell.effect)) 
    arrange(bcell.effect)
  jbins.bcell <- jsort.bcell$bin[1:keeptop]
  
  jsort.granu <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(desc(granu.effect))
    arrange(granu.effect)
  jbins.granu <- jsort.granu$bin[1:keeptop]
  
  jsort.eryth <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(desceryth.effect)) 
    arrange(eryth.effect)
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
} else {
  jsort.hspcs <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(mean.effect)
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
  
  jsort.bcell <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(desc(bcell.effect))
  jbins.bcell <- jsort.bcell$bin[1:keeptop]
  
  jsort.granu <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(desc(granu.effect))
  jbins.granu <- jsort.granu$bin[1:keeptop]
  
  jsort.eryth <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(desc(eryth.effect))
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
}



bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

bins.keep.lst <- list("Eryths" = jbins.eryth,
                      "Bcells" = jbins.bcell,
                      "Granulocytes" = jbins.granu,
                      "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames




# Plot effects: k9me3 on x axis  ------------------------------------------


jmark <- "H3K27me3"
ggplot(params.lst[[jmark]] %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) +
  geom_density() + 
  facet_wrap(~param) + 
  theme_bw() + 
  ggtitle(jmark) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# parms.keep

celltypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", celltypes.keep, ".Estimate", sep = "")

dat.params.wide.lst <- lapply(params.lst, function(jdat){
  jdat.tmp <- jdat %>% 
    filter(param %in% params.keep) %>% 
    group_by(bin) %>%
    filter(max(abs(estimate)) < 5)
  reshape2::dcast(data = jdat.tmp, formula = bin ~ param + mark, value.var = "estimate")
})
dat.params.wide.joined <- Reduce(f = left_join, x = dat.params.wide.lst[c("H3K4me1", "H3K4me3", "H3K27me3")], init = dat.params.wide.lst[["H3K9me3"]])


# jmarkref <- "H3K9me3"
# jmarkref <- "H3K4me1"


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/other_marks"
for (jmarkref in jmarks){
  outpdf <- file.path(outdir, paste0("H3K9me3_bins_by_pval.gene_gene_correlations.reference_", jmarkref, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
  
  mall <- ggplot(params.lst[[jmarkref]] %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) +
    geom_density() + 
    facet_wrap(~param) + 
    theme_bw() + 
    ggtitle(jmarkref) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(mall)
  
  
  jmarksfilt <- jmarks[jmarks != jmarkref]
  print(jmarksfilt)
  for (jparam in params.keep){
    mlst <- lapply(jmarksfilt, function(jmark){
      m <- ggplot(dat.params.wide.joined %>% filter(bin %in% k9.bins.names), 
                  aes_string(x = paste0(jparam, "_", jmarkref), y = paste0(jparam, "_", jmark))) + 
        geom_point() + 
        geom_density_2d() + 
        theme_bw() + 
        ggtitle(jmark, jparam) +
        geom_hline(yintercept = 0, linetype = "dotted") + 
        geom_vline(xintercept = 0, linetype = "dotted") + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      return(m)
    })
    JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], cols = 3)
  }
  dev.off()
}





