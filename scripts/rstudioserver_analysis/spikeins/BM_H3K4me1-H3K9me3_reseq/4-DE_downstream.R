# Jake Yeung
# Date of Creation: 2020-12-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/4-DE_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

dat.params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.rdata <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again/poisson_fit_bins.", jmark, ".2020-12-12.newannot2.witherrors.MoreBins.RData"))
  print(inf.rdata)
  load(inf.rdata, v=T)
  jfits.lst1 <- jfits.lst
  
  params.lst1 <- lapply(jfits.lst1, function(x){
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    x[xkeep]
  })

  # pvals.vec <- unlist(lapply(jfits.lst1, function(x) x$pval))
  pvals.dat <- data.frame(bin = names(jfits.lst1), pval = unlist(lapply(jfits.lst1, function(x) x$pval)), stringsAsFactors = FALSE)
  
  jnames <- names(jfits.lst1)
  names(jnames) <- jnames
  
  params.dat1 <- lapply(jnames, function(jname){
    jparams <- params.lst1[[jname]]
    if (is.null(jparams)){
      return(data.frame(NULL))
    }
    # assertthat::assert_that(!is.null(jparams))
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows() %>%
    left_join(., pvals.dat)
  if (jmark == "H3K9me3"){
   params.dat1 <- params.dat1 %>% 
    mutate(param = gsub("Eryth", "Eryths", param),
           param = gsub("Lymphoid", "Bcells", param))
  }
  return(params.dat1)
  
})

jmark <- "H3K4me1"
ggplot(dat.params.lst[[jmark]] %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) +
       geom_density() + 
       facet_wrap(~param) + 
       theme_bw() + 
       geom_vline(xintercept = 0, linetype = "dotted") + 
       theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.params.wide.lst <- lapply(dat.params.lst, function(jdat){
  reshape2::dcast(data = jdat, formula = bin ~ param + mark, value.var = "estimate")
})

dat.params.wide.joined <- Reduce(f = full_join, x = dat.params.wide.lst)

# fill NAs with zeros 
# dat.params.wide.joined[is.na(dat.params.wide.joined)] <- 0

ggplot(dat.params.wide.joined, aes(x = ClusterEryths.Estimate_H3K27me3, y = ClusterEryths.Estimate_H3K9me3)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  geom_density_2d() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.params.wide.joined, aes(x = ClusterBcells.Estimate_H3K4me1, y = ClusterBcells.Estimate_H3K27me3)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  geom_density_2d() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Compare log2FC across all celltypes but keep only signiicant pva --------

dat.merge1 <- left_join(dat.params.lst$H3K4me1, dat.params.lst$H3K27me3, by = c("bin", "param"))


ggplot(dat.merge1 %>% 
         filter(abs(estimate.x) < 5 & abs(estimate.y) < 5) %>%
         filter(pval.x < 1e-10), 
       aes(x = estimate.x, y = estimate.y)) + 
  geom_point(alpha = 0.25) + 
  geom_density_2d() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  facet_wrap(~param) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# # Check TSS bins used for matrix ------------------------------------------
# 
# inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")
# dat.genes <- fread(inf.genes)
# 



# 
# inf.rdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.H3K4me1.2020-12-11.newannot2.witherrors.RData")
# load(inf.rdata, v=T)
# 
# length(jfits.lst)
# 
# 
# 
# ggplot(params.dat1 %>% filter(abs(estimate1) < 5), aes(x = estimate1, fill = param)) +
#   geom_density() + 
#   facet_wrap(~param) + 
#   theme_bw() + 
#   geom_vline(xintercept = 0, linetype = "dotted") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


