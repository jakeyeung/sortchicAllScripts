# Jake Yeung
# Date of Creation: 2019-05-10
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/5-visualize_FACS_on_umap.R
# Visualize FACS on umap 

rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(hash)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxB6.R")



# Add cell counts to UMAP  ------------------------------------------------

# cellcounts
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1"); names(jmarks) <- jmarks

inf.csum <- "/Users/yeung/data/scchic/robjs/B6_objs/cell_sums_from_LDA_input.RData"
load(inf.csum, verbose = TRUE)

# umaps
inf.dats <- lapply(jmarks, function(jmark) paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.", jmark, ".RData"))
dat.umap.long <- lapply(inf.dats, LoadUmap) %>%
  bind_rows()

# merge
dat.merge <- left_join(dat.umap.long, count.sum)
dat.merge$cellsum.log10 <- log10(dat.merge$cellsum)

# jmarks <- "H3K4me1"
lapply(jmarks, function(jmark) PlotXYWithColor(dat.merge %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cellsum.log10", jtitle = jmark, jcol.low = "darkred", jcol.mid = "gray", jcol = "darkblue"))

