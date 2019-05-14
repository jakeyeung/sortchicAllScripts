# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/7b-create_clusters_from_traj_write_objects_H3K4me1.R
# From trajectory

rm(list=ls())

library(JFuncs)
library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)
library(ggrepel)

library(hash)
library(igraph)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Constants ---------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_from_traj"
dir.create(outdir)


# Load data ---------------------------------------------------------------

jmark <- "H3K4me1"
inf <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_", jmark, ".Rdata")
load(inf, v=T)

jtrajs <- names(trajs[[jmark]])
print(jtrajs)
# remove Tcell, it overlaps with lymphoid, use TcellIsland instead
jtrajs.filt <- jtrajs[which(jtrajs != "Tcell")]

for (jtraj in jtrajs.filt){
  dat.sub <- trajs[[jmark]][[jtraj]] %>% dplyr::select(cell)
  dat.sub$fname <- sapply(dat.sub$cell, function(x) paste0(x, ".sorted.bam"))
  outf.sub <- file.path(outdir, paste0(paste("bamlist", jtraj, jmark, sep = "-"), ".txt"))
  data.table::fwrite(dat.sub %>% dplyr::select(fname), file = outf.sub, sep = "\t", col.names = FALSE)
}
