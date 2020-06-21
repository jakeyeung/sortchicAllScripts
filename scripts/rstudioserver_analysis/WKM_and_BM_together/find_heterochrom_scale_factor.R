# Jake Yeung
# Date of Creation: 2020-06-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/find_heterochrom_scale_factor.R
# Find scale factor that helps normalize everything? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jctypes <- c("HSPCs", "Bcells", "Granulocytes", "Erythroblasts"); names(jctypes)

dat.stats <- lapply(jctypes, function(jctype){
  lapply(jmarks, function(jmark){
    indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts.countsByCluster"
    inf.hetero <- file.path(indir, paste0("hetero_and_totalcuts.", jctype, ".", jmark, ".txt"))
    print(inf.hetero)
    dat <- fread(file = inf.hetero)
    dat$mark <- jmark
    dat$ctype <- jctype
    return(dat)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

ggplot(subset(dat.stats), mapping = aes(x = cluster, fill = mark, y = ncuts.hetero / ncuts.total)) + geom_col() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, nrow = 1)

ggplot(subset(dat.stats), mapping = aes(x = cluster, fill = mark, y = ncell / ncuts.total)) + geom_col() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, nrow = 1)

ggplot(subset(dat.stats), mapping = aes(x = cluster, fill = mark, y = log10(ncell / ncuts.total))) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, nrow = 1)

ggplot(subset(dat.stats), mapping = aes(x = cluster, fill = mark, y = log10(ncuts.total / ncell))) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, nrow = 1)
