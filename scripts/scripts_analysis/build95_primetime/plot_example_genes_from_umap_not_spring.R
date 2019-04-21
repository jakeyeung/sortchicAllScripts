# Jake Yeung
# Date of Creation: 2019-04-19
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/plot_example_genes_from_umap_not_spring.R
# Plot example genes

rm(list=ls())


library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)


source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")


# Load --------------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"

load(inf, v=T)



# Constants ---------------------------------------------------------------


jmark <- "H3K4me1"

dat.umap <- dat.umap.long.new.lst[[jmark]]


jsize <- 2
jcols <- c("gray85", "gray50", "darkred")
jcols.gene <- jcols
jcols.gene[[3]] <- "darkblue"


# Check topisc ------------------------------------------------------------



top.dat <- data.frame(cell = rownames(tm.result.lst[[jmark]]$topics), as.data.frame(tm.result.lst[[jmark]]$topics))

top.dat <- left_join(top.dat, dat.umap)

jcnames <- paste("X", seq(30), sep = "")

mlst <- lapply(jcnames, function(jcname) PlotXYWithColor(top.dat, xvar = "umap1", yvar = "umap2", cname = jcname, jtitle = jcname))

pdf("/tmp/topics_H3K4me1.pdf", useDingbats = FALSE)
lapply(mlst, function(m) print(m))
dev.off()

# interesting topics 
jtops <- c(24, 23, 2)

# get genes

# Plot example genes ------------------------------------------------------

jhits <- c("Sox6", "Hbb-b2", "Irf8", "Irf4", "Foxo1", "S100a8", "Inpp4b", "Inhba", "Mef2c", "Car10", "Hs3st5", "Cd19", "Kit1", "Hoxa1", "Meis1", "Cd11b", "Cd3e")

# bcell markers

# jhits <- c("Cebpa", "Chi3l3", "Mac1", "Cd11b", "")
jhits <- c("Ets1", "Fli1", "Runx1", "Irf2", "Irf8", "Foxo1", "Ebf1", "Il7r", "Il2ra", "Prf1", "Hbb-b2", "Klf1", "Inpp4b", "S100a8", "Hs3st5", "Meis1", "Hoxa1", "Cd19")

jpeaks <- lapply(jhits, function(jhit){
  annots.filt <- GetPeaksFromGene(jhit, annots.lst[[jmark]], dist = 1000)
  jpeak <- SelectBestPeak(annots.filt$peaks, annots.lst[[jmark]], tm.result.lst[[jmark]])
  return(jpeak)
})

jcounts <- count.imputed.lst[[jmark]][unlist(jpeaks), ]
jcounts.dat <- data.frame(bin = rownames(jcounts), jcounts)
jcounts.dat <- tidyr::gather(jcounts.dat, key = cell, value = gene.counts, -bin)
dat.umap.with.counts <- left_join(jcounts.dat, dat.umap)

regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd),
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)

pdf("~/data/scchic/pdfs/H3K4me1_top_genes.pdf", useDingbats = FALSE)
for (i in seq(length(jpeaks))){
  
  jpeak <- jpeaks[[i]]
  jgene <- jhits[[i]]
  jsub <- subset(dat.umap.with.counts, bin == jpeak)
  jsub <- RankOrder(jsub, "gene.counts")
  # print(head(jsub))
  mdpt <- min(jsub$gene.counts) + (max(jsub$gene.counts) - min(jsub$gene.counts)) / 2
  mtop.gene <- ggplot(jsub, aes(x = umap1, y = -1 * umap2, color = gene.counts, order = orderrank)) + 
    geom_point(size = jsize) + 
    scale_color_gradient2(low = jcols.gene[[1]], mid = jcols.gene[[2]], high = jcols.gene[[3]], lim = c(min(jsub$gene.counts), max(jsub$gene.counts)), midpoint = mdpt) + 
    xlab("") + ylab("") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.ticks=element_blank(), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.border=element_blank()) + 
    ggtitle(jgene, jpeak)
  print(mtop.gene)
  
}
dev.off()
