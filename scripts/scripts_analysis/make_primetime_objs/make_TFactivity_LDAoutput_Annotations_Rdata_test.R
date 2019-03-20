# Jake Yeung
# Date of Creation: 2019-03-19
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_TFactivity_LDAoutput_Annotations_Rdata_test.R
# Test the 7 objects can do things


rm(list=ls())

library(ggplot2)
# library(ggrepel)

library(dplyr)
# library(hash)

library(JFuncs)
library(scales)
library(umap)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")


load("~/data/scchic/robjs/TFactivity_genelevels_objects.RData", v=T)


# Plot --------------------------------------------------------------------



jmarks.all <- list("H3K4me1" = "H3K4me1", "H3K4me3" = "H3K4me3", "H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")


jcolvec <- c("gray80", "gray50", "darkblue")
jmotif <- "Cebpb"
plts.lst <- lapply(jmarks.all, function(jmark){
  return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec))
})
multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)


names(count.mat.lst) <- jmarks.all

# plot a gene
ref.mark <- "H3K4me1"
jgene <- "Hoxc13"
out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.mark]])
(jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]]))
jscale.fac <- 1
jpseudo <- 10^-6
jsize <- 1

jmarks.sub <- c("H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")
system.time(
  PlotUmapAllMarks(jmarks.sub, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, custom.settings.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)
)



