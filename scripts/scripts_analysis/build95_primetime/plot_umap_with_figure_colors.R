# Jake Yeung
# Date of Creation: 2019-04-05
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/plot_umap_with_total_counts.R
# Plot umap with total counts

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

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")


# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"
load(inf, v=T)

# load trajectories
inf.trajs <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
load(inf.trajs, v=T)


# Constants ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks
jscales <- c(-1, -1, 1, 1)


jsize <- 2
colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K9me3 = "red1", H3K27me3 = "darkorange1")

jtops <- c(26, 13, 30, 23, 17, 24, 5, 7, 8, 12, 15, 18, 29, 28, 27)

jhits <- c("Sox6", "Hbb-b2", "Irf8", "Irf4", "Foxo1", "S100a8", "Inpp4b", "Inhba", "Mef2c", "Car10", "Hs3st5", "Cd19", "Kit1", "Hoxa1", "Meis1", "Cd11b", "Cd3e")

jcols <- c("gray85", "gray50", "darkred")
jcols.gene <- jcols
jcols.gene[[3]] <- "darkblue"




# Get counts total and append to UMAP  ------------------------------------




jmark <- jmarks[[1]]
jmark <- jmarks[[2]]
jmark <- jmarks[[3]]
jmark <- jmarks[[4]]
jscale <- 1


library(ggrastr)
umap.plots <- mapply(function(dat.umap, jcol, jscale){
  m <- ggplot(dat.umap, aes(x = umap1, y = jscale * umap2)) +
    # geom_point_rast(size = 0.5, alpha = 0.5, color = jcol) +
    geom_point(size = 0.5, alpha = 0.2, color = jcol) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.ticks=element_blank(), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.border=element_blank()) +
    # xlab("umap1") + ylab("umap2")
    xlab("") + ylab("")
  return(m)
}, dat.umap.long.new.lst, colsvec, jscales, SIMPLIFY = FALSE)

# plot interesting topics 
jmark <- "H3K4me1"
topics.mat <- tm.result.lst[[jmark]]$topics
dat.umap <- dat.umap.long.new.lst[[jmark]]

# add topic levels

tops.sub <- topics.mat[, jtops]
tops.long <- data.frame(cell = rownames(tops.sub), tops.sub) %>%
  tidyr::gather(key = "Topic", value = "Weight", -cell)
dat.umap.with.weights <- left_join(tops.long, dat.umap)


# plot top hits

jpeaks <- lapply(jhits, function(jhit){
  annots.filt <- GetPeaksFromGene(jhit, annots.lst[[jmark]], dist = 1000)
  jpeak <- SelectBestPeak(annots.filt$peaks, annots.lst[[jmark]], tm.result.lst[[jmark]])
  return(jpeak)
})

jcounts <- count.imputed.lst[[jmark]][unlist(jpeaks), ]
jcounts.dat <- data.frame(bin = rownames(jcounts), jcounts)
jcounts.dat <- tidyr::gather(jcounts.dat, key = cell, value = gene.counts, -bin)

dat.umap.with.counts <- left_join(jcounts.dat, dat.umap)



# plot trajectories

trajs.mark <- trajs[[jmark]]
ctypes <- c("eryth", "granu", "mega", "bcell", "tcell", "nkcell")
traj.jsize <- 1

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#C0C0C0")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")
m.trajs <- ggplot(dat.umap %>% mutate(umap2 = -1 * umap2), aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values=cbPalette)
for (ctype in ctypes){
  m.trajs <- m.trajs + geom_path(data = trajs.mark[[ctype]], color = "gray25", size = traj.jsize, alpha = 0.5)
}
print(m.trajs)


# Plot topic along trajectories 

# add topic counts to trajectories

# ctype <- ctypes[[1]]
# ctype <- ctypes[[1]]

pdf(paste0("~/data/scchic/pdfs/umap_with_colors_and_trajs_FIG2.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  multiplot(umap.plots[[1]], umap.plots[[2]], umap.plots[[3]], umap.plots[[4]], cols = 4)
  
  for (jtop.num in jtops){
    jtop <- paste0("X", jtop.num)
    jsub <- subset(dat.umap.with.weights, Topic == jtop)
    jsub <- RankOrder(jsub, "Weight")
    mdpt <- min(jsub$Weight) + (max(jsub$Weight) - min(jsub$Weight)) / 2
    mtop <- ggplot(jsub, aes(x = umap1, y = -1 * umap2, color = Weight, order = orderrank)) + geom_point(size = jsize) + 
      scale_color_gradient2(low = jcols[[1]], mid = jcols[[2]], high = jcols[[3]], lim = c(0, max(jsub$Weight)), midpoint = mdpt) + 
      xlab("") + ylab("") + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            panel.border=element_blank()) + 
      ggtitle(jtop)
    print(mtop)
  }
  
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
    mtop.gene.adj <- ggplot(jsub, aes(x = umap1, y = -1 * umap2, color = gene.counts, order = orderrank)) + 
      geom_point(size = jsize) + 
      scale_color_gradient2(low = jcols.gene[[1]], mid = jcols.gene[[2]], high = jcols.gene[[3]], lim = c(min(jsub$gene.counts), max(jsub$gene.counts)), midpoint = mdpt * 0.96) + 
      # scale_color_gradient2(low = jcols.gene[[1]], mid = jcols.gene[[2]], high = jcols.gene[[3]], lim = c(min(jsub$gene.counts), max(jsub$gene.counts)), midpoint = mdpt * 0.92) + 
      xlab("") + ylab("") + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            panel.border=element_blank()) + 
      ggtitle(jgene, jpeak)
    print(mtop.gene.adj)
  }
  print(m.trajs)
  # plot weights along trajectories
  for (ctype in ctypes){
    trajs.with.topweights <- left_join(tops.long, trajs.mark[[ctype]]) %>%
      filter(!is.na(lambda))
    m.traj.weights <- ggplot(trajs.with.topweights, aes(x = lambda, y = Weight)) + geom_point(alpha = 0.3) + 
      ggtitle(ctype) + facet_wrap(~Topic) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m.traj.weights)
  }
dev.off()
