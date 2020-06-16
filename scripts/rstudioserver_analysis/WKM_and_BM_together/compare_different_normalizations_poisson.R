# Jake Yeung
# Date of Creation: 2020-06-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/compare_different_normalizations_poisson.R
# Compare different normalizations for Poisson 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# Load poiisson outputs (wrangled) ----------------------------------------

jnorm <- "ByHetero"
inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm/fit_poisson_model_on_TSS.BM.NormMeth_", jnorm, ".RData")
load(inf, v=T)

jfits.dat.hetero.lst <- lapply(jmarks, function(jmark){
  jfit.dat <- do.call(rbind, jfits.lst.bymark[[jmark]])
  jfit.dat$bin <- rownames(jfit.dat)
  return(jfit.dat)
})

jfits.long.hetero.lst <- lapply(jmarks, function(jmark){
  jlong <- jfits.dat.hetero.lst[[jmark]] %>%
    dplyr::select(-c(dev.diff, df.diff)) %>%
    reshape2::melt(., id.vars = c("bin", "pval"), variable.name = "cluster", value.name = "logLambda")
  jlong$mark <- jmark
  # extract the X intercept as separate column 
  jlong.noint <- subset(jlong, cluster != "X.Intercept.")
  jlong.int <- subset(jlong, cluster == "X.Intercept.")  %>%
    dplyr::rename(logintercept = logLambda) %>%
    dplyr::select(bin, logintercept, mark)
  jlong.merge <- left_join(jlong.noint, jlong.int)
  return(jlong.merge)
})


jnorm2 <- "ByTotalFromBins"
inf2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm/fit_poisson_model_on_TSS.BM.NormMeth_", jnorm2, ".RData")
# inf2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression/fit_poisson_model_on_TSS.2020-06-05.RData"
load(inf2, v=T)


jfits.dat.total.lst <- lapply(jmarks, function(jmark){
  jfit.dat <- do.call(rbind, jfits.lst.bymark[[jmark]])
  jfit.dat$bin <- rownames(jfit.dat)
  return(jfit.dat)
})

jfits.long.total.lst <- lapply(jmarks, function(jmark){
  jlong <- jfits.dat.total.lst[[jmark]] %>%
    dplyr::select(-c(dev.diff, df.diff)) %>%
    reshape2::melt(., id.vars = c("bin", "pval"), variable.name = "cluster", value.name = "logLambda")
  jlong$mark <- jmark
  # extract the X intercept as separate column 
  jlong.noint <- subset(jlong, cluster != "X.Intercept.")
  jlong.int <- subset(jlong, cluster == "X.Intercept.")  %>%
    dplyr::rename(logintercept = logLambda) %>%
    dplyr::select(bin, logintercept, mark)
  jlong.merge <- left_join(jlong.noint, jlong.int)
  return(jlong.merge)
})

# plot outputs of two normalizations --------------------------------------

jparam <- "ClusterBcells"
# jparam <- "X.Intercept."
jmark <- "H3K4me3"

bins.filt <- subset(jfits.dat.total.lst[[jmark]], X.Intercept. > -15)$bin

jmerge <- full_join(subset(jfits.dat.total.lst[[jmark]], select = c(bin, ClusterBcells, ClusterErythroblasts, ClusterGranulocytes, X.Intercept.)), 
                    subset(jfits.dat.hetero.lst[[jmark]], select = c(bin, ClusterBcells, ClusterErythroblasts, ClusterGranulocytes, X.Intercept.)), by = "bin")

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterBcells.x, y = ClusterBcells.y))  + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterBcells.x - ClusterBcells.y))  + 
  geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterGranulocytes.x - ClusterGranulocytes.y))  + 
  geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterErythroblasts.x - ClusterErythroblasts.y))  + 
  geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterBcells.x)) + 
  geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterBcells.y)) + 
  geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterGranulocytes.x)) + 
  geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterGranulocytes.y)) + 
  geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterErythroblasts.x)) + 
  geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)

ggplot(jmerge %>% filter(bin %in% bins.filt & ClusterBcells.x > -10), aes(x = ClusterErythroblasts.y)) + 
  geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)

# plot(density(jfits.dat.total.lst$))

# check pax5 H3K27me3 
jmark <- "H3K27me3"
jgene <- "Hbb-y"
jgene <- "Tal1"
jgene <- "Hbb-bh2"

(jbin <- subset(jfits.dat.total.lst[[jmark]], grepl(jgene, bin))$bin[[1]])


head(subset(jfits.dat.total.lst[[jmark]], bin == jbin))
head(subset(jfits.dat.hetero.lst[[jmark]], bin == jbin))


# Plot density?  ----------------------------------------------------------

plot(density(jfits.dat.total.lst$H3K4me1[[jparam]])); abline(v = 0)
plot(density(jfits.dat.total.lst$H3K4me3[[jparam]])); abline(v = 0)
plot(density(jfits.dat.total.lst$H3K27me3[[jparam]])); abline(v = 0)
# plot(density(jfits.dat.total.lst[[jmark]][[jparam]]))

