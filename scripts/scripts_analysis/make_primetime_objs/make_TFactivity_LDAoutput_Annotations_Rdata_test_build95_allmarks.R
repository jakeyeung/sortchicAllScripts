# Jake Yeung
# Date of Creation: 2019-03-26
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_TFactivity_LDAoutput_Annotations_Rdata_test_build95.R
# description

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


load("~/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks.RData", v=T)


# Plot --------------------------------------------------------------------



jmarks.all <- c("H3K4me1" = "H3K4me1", "H3K4me3" = "H3K4me3", "H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")

jmarks.sub <- jmarks.all

jcolvec <- c("gray80", "gray50", "darkblue")
jmotif <- "Tal1"
jmotif <- "Tfdp1"
jmotif <- "Pml"



jmotif <- "Bptf"
jmotif <- "Bcl3"
jmotif <- "Irf1"

jmotif <- "Sox6"



jmotif <- "Hoxc6"
jmotif <- "Cebpb"

jmotif <- "Stat6"
plts.lst <- lapply(jmarks.all, function(jmark){
  return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 1.5, colvec = jcolvec))
})
multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)

# names(count.mat.lst) <- jmarks.all


# plot a gene
ref.mark <- "H3K27me3"
jgene <- "Hoxc6"
jgene <- "Stat6"

out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.mark]])
(jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]]))
jscale.fac <- 1
jpseudo <- 10^-6
jsize <- 1

jmarks.sub <- c("H3K27me3" = "H3K27me3", "H3K9me3" = "H3K9me3")
system.time(
  PlotUmapAllMarks(jmarks.sub, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, custom.settings.new.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)
)

# plot zscores?
head(mara.outs$H3K27me3$zscores)
head(mara.outs$H3K27me3$zscores)

zscores.merge <- left_join(mara.outs$H3K27me3$zscores, mara.outs$H3K9me3$zscores, by = "motif") %>%
  mutate(motif.lab = ifelse(zscore.x >= 0.6 | zscore.y >= 6, motif, NA))

ggplot(zscores.merge, aes(x = zscore.x, y = zscore.y, label = motif.lab)) + ggrepel::geom_text_repel(size = 6) + 
  geom_point() + theme_bw() + xlab("H3K9me3 Zscore") + ylab("H3K27me3 Zscore")
