# Jake Yeung
# Date of Creation: 2020-06-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/check_lm_var_correction_output.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# Load gene sets ----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

# Load inputs -------------------------------------------------------------

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits_no_var.2020-06-26.RData"
load(inf, v=T)

dat.params.all.lst <- lapply(dat.params.all.lst, function(jdat){
  jsub.int <- subset(jdat, params.mean == "clstHSC")
  jsub.int.tomerge <- data.frame(rname = jsub.int$rname, intercept = jsub.int$jmean)
  jdat.merged <- left_join(jdat, jsub.int.tomerge)
})

# Check DE of gene sets  --------------------------------------------------

dat.params.all.annot <- lapply(jmarks, function(jmark){
  jdat <- dat.params.all.lst[[jmark]]
  jdat$gene <- sapply(jdat$rname, function(x) strsplit(x, ";")[[1]][[2]])
  jdat$ens <- sapply(jdat$gene, AssignHash, g2e.hash)
  jdat$mark <- jmark
  jdat <- subset(jdat, !is.na(ens))
  return(jdat)
}) %>%
  bind_rows()


# do by gene lists 
jens.vec.names <- names(de.ens.sorted.stringent)
names(jens.vec.names) <- jens.vec.names

params.bygsets <- lapply(jens.vec.names, function(jname){
  ens.vec <- de.ens.sorted.stringent[[jname]]
  jsub <- subset(dat.params.all.annot, ens %in% ens.vec)
  jsub$gset <- jname
  jsub <- jsub %>%
    filter(params.mean != "clstHSC") %>%
    mutate(logfc = jmean - intercept)
  return(jsub)
}) %>%
  bind_rows()

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(params.bygsets, aes(x = logfc, fill = gset)) + facet_grid(mark~params.mean) + geom_density(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) + scale_fill_manual(values = cbPalette)

for (jmark in jmarks){
  m <- ggplot(params.bygsets %>% filter(mark == jmark), aes(x = logfc, fill = gset)) + facet_grid(gset~params.mean) + geom_density(alpha = 0.25) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_vline(xintercept = 0) + scale_fill_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}

# show genomewide

m <- ggplot(dat.params.all.annot %>% filter(params.mean != "clstHSC"), aes(x = jmean - intercept, fill = mark)) + facet_wrap(~params.mean) + geom_density(alpha = 0.25) + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0)
print(m)


