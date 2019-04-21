# Jake Yeung
# Date of Creation: 2019-04-21
# File: ~/projects/scchic/scripts/scripts_analysis/tf_activity_debugging/explore_LDA_100kb_bins_near_genes.R
# In peak analysis, nearby peaks can have very different spatial distributions. Is this true for 100kb bins?

rm(list=ls())

library(dplyr)
library(ggplot2)


source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/Aux.R")

# Load  -------------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"

load(inf, v=T)

bsize <- 100000


# Get all peaks near a gene -----------------------------------------------

jmark <- "H3K4me1"
jgene <- "Cebpb"
cvec <- c("gray85", "gray50", scales::muted("darkblue"))

annots.sub <- annots.lst[[jmark]] %>% filter(SYMBOL == jgene) %>% filter(distanceToTSS == min(abs(distanceToTSS)))

jchromo <- annots.sub$seqnames[[1]]
jchromo.str <- paste0("^", jchromo, ":")

counts.imputed <- t(tm.result.lst[[jmark]]$topics %*% tm.result.lst[[jmark]]$terms)
counts.mat.sub <- counts.imputed[grepl(jchromo.str, rownames(counts.imputed)), ]
counts.mat.sub.long <- tidyr::gather(data = as.data.frame(counts.mat.sub) %>% mutate(bin = rownames(counts.mat.sub)) %>% rowwise() %>% mutate(start = as.numeric(GetStart(bin))),
                                     key = "cell", value = "exprs", c(-bin, -start)) %>%
  arrange(start)
# add umap
counts.mat.sub.long <- left_join(counts.mat.sub.long, dat.umap.long.trajs[[jmark]] %>% dplyr::select(umap1, umap2, cell, louvain))

# test
# annotpeak <- "chr4:115020000-115120000"
annotpeak <- annots.sub$region_coord[[1]]
PlotImputedPeaks3(counts.mat.sub.long %>% filter(bin == annotpeak), peaks.keep = annotpeak, jmark, gname = jgene, jsize = 1, jcolvec = cvec, .log = FALSE)

# do all near tal1
jpeaks <- unique(counts.mat.sub.long$bin)

jpeaks.i <- seq(which(jpeaks == annotpeak) - 50, which(jpeaks == annotpeak) + 50)

jpeaks.filt <- jpeaks[jpeaks.i]

jgenemid <- as.numeric(GetStart(annotpeak)) + bsize / 2
pdf(paste0("/tmp/plot_bins_near_", jgene, ".pdf"), useDingbats = FALSE)
for (jpeak in jpeaks.filt){
  print(jpeak)
  peakmid <- as.numeric(GetStart(jpeak)) + bsize / 2
  jdist <- peakmid - jgenemid
  print(PlotImputedPeaks3(counts.mat.sub.long, jpeak, jmark, gname = paste(jgene, jdist), jsize = 1, jcolvec = cvec, .log = FALSE)) 
}
dev.off()

