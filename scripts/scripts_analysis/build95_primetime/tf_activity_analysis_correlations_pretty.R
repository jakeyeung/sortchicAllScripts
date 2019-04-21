# Jake Yeung
# Date of Creation: 2019-04-15
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/tf_activity_analysis_correlations_pretty.R
# Make pretty 

rm(list=ls())

tstart <- Sys.time()

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(JFuncs)
library(umap)
library(ggrepel)
library(tidyr)

source("scripts/Rfunctions/MaraDownstream.R")

source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"
inf.pvals <- "/Users/yeung/data/scchic/robjs/TFactivity_zscore_analysis.RData"
# inf.trajs <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
# inf spring umaps 
inf.spring <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"


# Functions ---------------------------------------------------------------

RenameTraj <- function(x){
  if (x == "eryth"){
    return("Erythryoid")
  } else if (x == "granu") {
    return("Myeloid")
  } else if (x == "lymphoid"){
    return("Lymphoid") 
  }
}

RenameTraj.rev <- function(x){
  if (x == "Erythryoid"){
    return("eryth")
  } else if (x == "Myeloid") {
    return("granu")
  } else if (x == "Lymphoid"){
    return("lymphoid") 
  }
}

# Constants ---------------------------------------------------------------

pcutoff <- -7.5

jscale.fac <- 10^6
jpseudo <- 0
jsize <- 1
make.plots <- FALSE
# make.plots <- TRUE


# Load data  --------------------------------------------------------------

load(inf, v=T)
load(inf.pvals, v=T)
load(inf.trajs, v=T)
load(inf.spring, v=T)

rm(count.mat.lst)

head(mara.outs$H3K4me1$zscore)




# Show top hits -----------------------------------------------------------

pvals.top <- pvals.long %>%
  group_by(mark) %>%
  dplyr::top_n(n = 10, wt = -log10pval) %>%
  arrange(log10pval)
print(dim(pvals.top))
print(split(pvals.top, f = pvals.top$mark))



# Calculate correlations --------------------------------------------------

ctypes <- c("granu", "lymphoid", "eryth")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

# ref.mark <- "H3K4me3"
ref.mark <- "H3K4me1"

# for (ref.mark in jmarks){
  jmark <- ref.mark
  
  dat.trajs.long <- dat.trajs.long %>% 
    filter(mark == jmark)
  
  # add x1 and x2 to umap1 and umap2 
  
  dat.umap.act.filt <- left_join(dat.merged.lst[[jmark]] %>% dplyr::select(-umap1, -umap2), dat.trajs.long %>% dplyr::select(X1, X2, mark, cell, louvain)) %>%
    dplyr::rename(umap1 = X1, umap2 = X2)
  
  # flip y axis optionally 
  umap.fac <- ifelse(ref.mark %in% c("H3K4me1", "H3K4me3"), 1, 1)
  
  if (make.plots){
    pdfout <- paste0("~/data/scchic/pdfs/TF_analysis_gene_correlations.WithSpringTrajs.refmark.", ref.mark, ".", Sys.Date(), ".Pretty.pdf")
    # assertthat::assert_that(!file.exists(pdfout))
    pdf(pdfout, useDingbats = FALSE)
  }
  
  
  if (nrow(subset(pvals.long, mark == ref.mark)) < 1){
    stop(paste("Ref markd oes not match anything:", ref.mark))
  }
  
  head(print(data.frame(subset(pvals.long, mark == ref.mark) %>% arrange(log10pval))), n = 50)
  
  jcolvec <- c("gray85", "gray50", "darkblue")
  jcolvec.motif <- c("gray85", "gray50", "red")
  
  peaks.filt <- (annots.lst[[jmark]] %>% filter(distanceToTSS == 0))$region_coord
  rows.i <- rownames(count.imputed.lst[[jmark]]) %in% peaks.filt
  counts.mat.sub <- count.imputed.lst[[jmark]][rows.i, ]
  
  
  
  
  
  
  
  
  # Do large-scale correlations in a trajectory-specific way  ---------------
  
  annots.sub <- annots.lst[[jmark]] %>% 
    filter(region_coord %in% rownames(counts.mat.sub)) %>%
    dplyr::select(SYMBOL, region_coord, distanceToTSS) %>%
    dplyr::rename(bin = region_coord,
                  motif = SYMBOL) %>%
    rowwise() %>%
    dplyr::mutate(start = as.numeric(GetStart(bin)),
                  end = as.numeric(GetEnd(bin)),
                  binmid = start + (end - start) / 2,
                  dist = distanceToTSS) %>%
    # pick best bin based on smallest distance
    group_by(motif) %>%
    dplyr::filter(dist == min(abs(dist)))
  
  # pick closest 
  counts.mat.sub.long <- tidyr::gather(data = as.data.frame(counts.mat.sub) %>% mutate(bin = rownames(counts.mat.sub)), 
                                       key = "cell", value = "exprs", -bin)
  # add trajectory info
  
  traj.annot <- lapply(ctypes, function(jtraj) data.frame(cell = trajs.spring[[jmark]][[jtraj]]$cell, traj = jtraj, lambda = trajs.spring[[jmark]][[jtraj]]$lambda)) %>%
    bind_rows()
  
  traj.annot <- left_join(traj.annot, dat.trajs.long %>% dplyr::select(X1, X2, cell)) %>%
    dplyr::rename(umap1 = X1, umap2 = X2)
  
  # assign bin to gene
  counts.mat.sub.long <- left_join(counts.mat.sub.long, annots.sub)
  
  # add trajectory information
  counts.mat.sub.long <- left_join(counts.mat.sub.long, traj.annot) %>%
    filter(!is.na(traj))
  
  # add activity info
  counts.mat.sub.long.inner <- inner_join(counts.mat.sub.long, dat.umap.act.filt %>% dplyr::select(-umap1, -umap2)) %>%
    group_by(motif) %>%
    mutate(exprs = scale(exprs, center = TRUE, scale = TRUE),
           activity = scale(activity, center = TRUE, scale = TRUE))
  
  
  # do genome-wide correlations? 
  counts.fit <- counts.mat.sub.long.inner %>%
    group_by(motif, traj) %>% 
    summarise(cor.out = cor(x = exprs, y = activity)) %>%
    arrange(desc(abs(cor.out)))
  
  # # do genome-wide linear model fits
  # counts.lmfit <- counts.mat.sub.long.inner %>%
  #   group_by(motif, traj) %>%
  #   summarise(lm.fit = GetSlope(x = exprs, y = activity)) %>%
  #   arrange(desc(abs(cor.out)))
  
  # add zscore and zscore pvalue information
  zscore.sub <- subset(zscore.long, mark == jmark & seed == 1) %>%
    dplyr::select(motif, mark, zscore.real)
  pval.sub <- subset(pvals.long, mark == jmark)
  
  counts.fit <- left_join(counts.fit, zscore.sub)
  counts.fit <- left_join(counts.fit, pval.sub)
  counts.fit %>% arrange(log10pval)
  
  counts.fit.filt <- counts.fit %>% 
    group_by(motif) %>% 
    filter(abs(cor.out) == max(abs(cor.out))) %>% 
    mutate(motif.lab = ifelse(log10pval < pcutoff, motif, NA))
  
  m.logpval.cor <- ggplot(counts.fit.filt, 
                          aes(y = -log10pval, x = cor.out, label = motif)) + geom_point() +  
    geom_text() + 
    facet_wrap(~traj) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Pearson Correlation") + 
    ylab("-log10(Pvalue of Zscore)") + 
    xlim(c(-1, 1))
  print(m.logpval.cor)
  
  m.zscore.cor <- ggplot(counts.fit.filt, 
                         aes(y = zscore.real, x = cor.out, label = motif)) + geom_point() +  
    geom_text() + 
    facet_wrap(~traj) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    xlim(c(-1, 1))  +
    xlab("Pearson Correlation") + 
    ylab("Motif Zscore")
  print(m.zscore.cor)
  
  # plot prettier with a few names and darkness by -log10 pvalue 
  m.logpval.cor <- ggplot(counts.fit.filt, 
                          aes(y = -log10pval, x = cor.out, label = motif.lab, color = zscore.real)) + 
    geom_point() +  
    geom_text_repel(alpha = 1) + 
    facet_wrap(~traj) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    xlim(c(-1, 1)) + 
    xlab("Pearson Correlation") + 
    ylab("-log10(Pvalue of Zscore)")  + 
    scale_color_gradient2(low = "gray85", mid = "gray50", high = scales::muted("darkblue"), midpoint = 1, lim = c(0, 2.75), name = "Zscore")
  print(m.logpval.cor)
  
  jmotifs <- (counts.fit.filt %>% filter(!is.na(motif.lab)) %>% arrange(log10pval))$motif
  
  # plot hits
  # jmotif <- "Sox30"
  # jmotif <- "Bcl3"
  
  for (jmotif in jmotifs){
    
    print(paste("Printing:", jmotif))
    
    jgene <- jmotif
    jsub <- subset(counts.mat.sub.long.inner, motif == jmotif)
    # rename trajectories
    

    jsub$traj <- sapply(jsub$traj, RenameTraj)
    ctypes.renamed <- unique(jsub$traj)
    (n.cells <- unique(length(jsub$cell)))
    
    m.exprsVactivity <- ggplot(jsub, aes(x = exprs, y = activity, color = traj)) + geom_point() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      xlab("scChIC Signal") + ylab("Motif Activity") + 
      ggtitle(paste(jmark, jmotif))
    
    m.exprsVactivity.split <- ggplot(jsub, aes(x = exprs, y = activity, color = traj)) + 
      geom_point(alpha = 0.5) + facet_wrap(~traj) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
      xlab("scChIC Signal") + ylab("Motif Activity") + 
      ggtitle(paste(jmark, jmotif))
    
    m.exprsVactivity.split.black <- ggplot(jsub, aes(x = exprs, y = activity)) + 
      geom_point(alpha = 0.5) + facet_wrap(~traj) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
      xlab("scChIC Signal") + ylab("Motif Activity") + 
      ggtitle(paste(jmark, jmotif))
    
    # plot each traj individually
    m.exprsVactivity.split.lst <- lapply(ctypes.renamed, function(ctype){
      pcor <- signif(subset(counts.fit, motif == jmotif & traj == RenameTraj.rev(ctype))$cor.out, 2)
      m.tmp <- ggplot(jsub %>% filter(traj == ctype), aes(x = exprs, y = activity)) + 
        geom_point(alpha = 0.5, color = "black") + facet_wrap(~traj) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_smooth(method = "lm", se = FALSE, color = "gray50") + 
        xlab("scChIC Signal") + ylab("Motif Activity") + 
        ggtitle(paste(jmark, jmotif, "Corr:", pcor))
    })
    
    m.pseudo <- ggplot(jsub %>% gather(key = "type", value = "exprsORact", c(exprs, activity)), aes(x = lambda, y = exprsORact, color = type)) + 
      geom_point(alpha = 0.5) + facet_wrap(~traj) + 
      ggtitle(paste(jmark, jmotif)) +
      xlab("Pseudotime") + ylab("Scaled scChICor Activity") + 
      scale_color_manual(labels = c("Motif Activity", "scChIC Signal"), values = c("darkred", "darkblue"), name = "") + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
    m.pseudo.lst <- lapply(ctypes.renamed, function(ctype){
      m.tmp <- ggplot(jsub %>% filter(traj == ctype) %>% gather(key = "type", value = "exprsORact", c(exprs, activity)), aes(x = lambda, y = exprsORact, color = type)) + 
        geom_point(alpha = 0.5) + facet_wrap(~traj) + 
        ggtitle(paste(jmark, jmotif)) +
        xlab("Pseudotime") + ylab("Scaled scChICor Activity") + 
        scale_color_manual(labels = c("Motif Activity", "scChIC Signal"), values = c("darkred", "darkblue"), name = "") + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    })
    
    print(m.exprsVactivity)
    print(m.exprsVactivity.split)
    print(m.exprsVactivity.split.black)
    print(m.exprsVactivity.split.lst)
    print(m.pseudo)
    print(m.pseudo.lst)
    
    
    # plot on UMAP 
    
    # plot positive or negative correlations
    out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.mark]], dist = 0)
    (jpeak <- SelectBestPeak(out.sub$peaks, NULL, tm.result.lst[[ref.mark]])[[1]])
    # get gene form bin
    # (jgene <- subset(annots.lst[[jmark]], region_coord == jbin)$SYMBOL)
    # m.peak <- PlotImputedPeaks2(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],  
    #                             # use.count.mat = count.mat.lst[[jmark]],
    #                             use.count.mat = NULL,
    #                             usettings=custom.settings.new.lst[[jmark]], 
    #                             gname = jgene,
    #                             jsize = jsize, jcolvec = jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo, y.axis.factor = -1)
    m.peak <- PlotImputedPeaks3(counts.mat.sub.long, 
                                jpeak, jmark, gname = jgene, jcolvec = jcolvec, .log = TRUE, jpseudo = 0, jscale = 10^7, jsize = jsize)
    m.motif <- PlotMotifInUmap(jmotif, dat.umap.act.filt %>% mutate(umap2 = umap.fac * umap2), mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec.motif)
    multiplot(m.peak, m.motif, cols = 2)
    
  }
  
  if (make.plots){
    dev.off()
  }
# }

print(Sys.time() - tstart)

# Below for testing -------------------------------------------------------


# 

# 
# jmotif <- "Stat6"
# jgene <- "Stat6"
# 
# jmotif <- "Ebf1"
# jgene <- "Ebf1"
# 
# jmotif <- "Foxc1"
# jgene <- "Foxc1"
# 
# jmotif <- "Pou2f2"
# jgene <- "Pou2f2"
# 
# jmotif <- "Gfi1"
# jgene <- "Gfi1"

# plot positive or negative correlations

# 
# # jmotif <- "Ebf1"
# jmotif <- "Cebpb"
# # jmotif <- "Foxc1"
# # 
# jgene <- jmotif
# # 
# # 
# out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.mark]], dist = 0)
# (jpeak <- SelectBestPeak(out.sub$peaks, NULL, tm.result.lst[[ref.mark]])[[1]])
# # get gene form bin
# # (jgene <- subset(annots.lst[[jmark]], region_coord == jbin)$SYMBOL)
# # m.peak <- PlotImputedPeaks2(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
# #                             # use.count.mat = count.mat.lst[[jmark]],
# #                             use.count.mat = NULL,
# #                             usettings=custom.settings.new.lst[[jmark]],
# #                             gname = jgene,
# #                             jsize = jsize, jcolvec = jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo, y.axis.factor = -1)
# m.peak <- PlotImputedPeaks3(counts.mat.sub.long,
#                             jpeak, jmark, gname = jgene, jcolvec = jcolvec, .log = TRUE, jpseudo = 0, jscale = 10^7, jsize = jsize)
# m.motif <- PlotMotifInUmap(jmotif, dat.umap.act.filt %>% mutate(umap2 = umap.fac * umap2), mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec)
# multiplot(m.peak, m.motif, cols = 2)
# 
# # add UMAP1 and UMAP2 to counts
# 
# m.peak <- PlotImputedPeaks3(counts.mat.sub.long,
#                              jpeak, jmark, gname = jgene, jcolvec = jcolvec, .log = TRUE, jpseudo = 0, jscale = 10^7, jsize = 0.5)
# print(m.peak)

# add 
# 
# 

# # what's correlation along the trajectories?
# 
# act.sub <- dat.umap.act.filt %>% filter(motif == jmotif) %>% select(cell, motif, activity)
# exprs.sub <- count.imputed.lst[[jmark]][jpeak, ]
# exprs.sub.dat <- data.frame(exprs = exprs.sub, cell = names(exprs.sub))
# merged.sub <- left_join(act.sub, exprs.sub.dat)
# 
# merged.sub.long <- tidyr::gather(data = merged.sub, key = "type", value = "exprsORactivity", -cell, -motif) %>%
#   group_by(type) %>%
#   mutate(exprsORactivity = scale(exprsORactivity, center = TRUE, scale = TRUE))

# # add trajectory 
# 
# for (ctype in ctypes){
#   act.exprs.traj <- left_join(trajs.spring[[jmark]][[ctype]], merged.sub.long)
#   # plot correlation
#   jtitle <- paste("Gene:", jgene, "Traj:", ctype)
#   m.cor <- ggplot(act.exprs.traj %>% spread(key = "type", value = "exprsORactivity"), aes(x = exprs, y = activity)) + geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_smooth(method = "lm", se = FALSE) + 
#     ggtitle(jtitle)
#   m <- ggplot(act.exprs.traj, aes(x = lambda, y = exprsORactivity, color = type)) + geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     ggtitle(jtitle) + geom_smooth(method = "lm", se = FALSE) + 
#     xlab("Pseudotime") + ylab("Scaled scChICor Activity")
#   print(m.cor)
#   print(m)
# }
# 
# 
# jmotif <- "Foxc1"
# jmotif <- "Cebpb"
# jmotif <- "Ebf1"
# jmotif <- "Foxc1"
# jmotif <- "Tal1"
# jsub <- subset(counts.mat.sub.long.inner, motif == jmotif)
# (n.cells <- unique(length(jsub$cell)))
# 
# ggplot(jsub, aes(x = exprs, y = activity, color = traj)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(paste(jmark, jmotif))
# 
# ggplot(jsub, aes(x = exprs, y = activity, color = traj)) + 
#   geom_point(alpha = 0.5) + facet_wrap(~traj) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_smooth(method = "lm", se = FALSE) + 
#   ggtitle(paste(jmark, jmotif))
# 
# # along trajecotry?
# ggplot(jsub %>% gather(key = "type", value = "exprsORact", c(exprs, activity)), aes(x = lambda, y = exprsORact, color = type)) + 
#   geom_point(alpha = 0.5) + facet_wrap(~traj) + 
#   ggtitle(paste(jmark, jmotif)) + 
#   xlab("Pseudotime") + ylab("Scaled scChICor Activity")
# theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
