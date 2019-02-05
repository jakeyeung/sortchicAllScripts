# Jake Yeung
# Date of Creation: 2019-01-30
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/debug.R
# Debug why UMAPs are different

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(forcats)
library(ggrepel)
library(biomaRt)

library(igraph)  # louvain

library(Gviz)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

source("scripts/Rfunctions/PlotFunctions.R")


load("/private/tmp/.RData", v=T)

plot(dat.umap$layout[, 1], dat.umap$layout[, 2])

dat.umap.osx <- umap(topics.mat, config = custom.settings)

par(pty="s")
plot(dat.umap.osx$layout[, 1], dat.umap.osx$layout[, 2], asp=1)

PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)
