# Jake Yeung
# Date of Creation: 2019-03-11
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/mara_downstream_both_marks_prom_enh_1kb_10kb.R
# Do 3 arrow analysis 

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

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")

# Functions ---------------------------------------------------------------


# Get dirs ----------------------------------------------------------------

# plotf <- "~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/mara/H3K4me1_and_H3K4me3_motifs_promenh.pdf"

jchips <- c("H3K4me1", "H3K4me3")

# types <- "all"
types <- c("enhancer", "all", "promoter")
zscore.thres <- 0.75

zscores.merged.lst <- lapply(types, function(type){
  
  mdirs <- lapply(jchips, function(jchip){
    if (type == "all"){
      mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                     jchip, 
                     ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                     "/",
                     "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jchip, ".filt_0.99.center_TRUE")
    } else if (type == "promoter"){
      mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                     jchip, 
                     ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.1000-",
                     "/",
                     "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jchip, ".filt_0.99.center_TRUE")
    } else if (type == "enhancer"){
      mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                     jchip, 
                     ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.10000-",
                     "/",
                     "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jchip, ".filt_0.99.center_TRUE")
    }
    assertthat::assert_that(dir.exists(mdir))
    return(mdir)
  })
  mara.outs <- lapply(mdirs, LoadMARA)
  
  act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long)
  zscores.merged <- left_join(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, by = "motif")
  zscores.merged.mean <- zscores.merged %>%
    mutate(zscore.mean = (zscore.x + zscore.y) / 2)
  
  
  zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
    ifelse(max(jrow[[2]], jrow[[3]]) > zscore.thres, jrow[[1]], NA)
  })
  zscores.merged$type <- type
  return(zscores.merged)
})

zscores.merged <- bind_rows(zscores.merged.lst)

# integrate the two datasets together and then show the arrows 

zscores.diff <- zscores.merged %>%
  group_by(motif) %>%
  summarise(zscore.x.enh = zscore.x[[1]], zscore.x.all = zscore.x[[2]], 
            zscore.y.enh = zscore.y[[1]], zscore.y.all = zscore.y[[2]],
            zscore.x.prom = zscore.x[[3]], zscore.x.prom = zscore.x[[3]],
            zscore.y.prom = zscore.y[[3]], zscore.y.prom = zscore.y[[3]],
              dist = sqrt((zscore.x.enh - zscore.x.all)^2 + (zscore.y.enh - zscore.y.all) ^ 2))

m.zscore <- ggplot(zscores.merged, aes(x = zscore.x, y = zscore.y, label = motif.lab)) + 
  geom_point(alpha = 0.5) + theme_bw() + # geom_text_repel() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab(paste("Motif Zscore:", jchips[[1]])) + ylab(paste("Motif Zscore:", jchips[[2]])) + 
  facet_wrap(~type) + 
  xlim(c(0, 2.5)) + ylim(c(0, 2.5))
print(m.zscore)

m.zscore.arrows <- ggplot() + 
  geom_segment(mapping = aes(x = zscore.x.all, xend = zscore.x.prom, y = zscore.y.all, yend = zscore.y.prom), 
               arrow = arrow(), data = zscores.diff %>% filter(zscore.x.all > 1 | zscore.y.all > 1)) + 
  geom_point(data = zscores.diff, aes(x = zscore.x.all, y = zscore.y.all)) + 
  geom_text_repel(data = zscores.diff %>% filter(zscore.x.all > 1 | zscore.y.all > 1), aes(label = motif, x = zscore.x.all, y = zscore.y.all)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Motif from all TFBS to removing 1kb from promoters")
print(m.zscore.arrows)

m.zscore.arrows.10kb <- ggplot() + 
  geom_segment(mapping = aes(x = zscore.x.prom, xend = zscore.x.enh, y = zscore.y.prom, yend = zscore.y.enh), 
               arrow = arrow(), data = zscores.diff %>% filter(zscore.x.all > 1 | zscore.y.all > 1)) + 
  geom_segment(mapping = aes(x = zscore.x.all, xend = zscore.x.prom, y = zscore.y.all, yend = zscore.y.prom), 
               data = zscores.diff %>% filter(zscore.x.all > 1 | zscore.y.all > 1)) + 
  geom_point(data = zscores.diff, aes(x = zscore.x.all, y = zscore.y.all)) + 
  geom_point(data = zscores.diff, aes(x = zscore.x.prom, y = zscore.y.prom)) + 
  geom_text_repel(data = zscores.diff %>% filter(zscore.x.all > 1 | zscore.y.all > 1), aes(label = motif, x = zscore.x.all, y = zscore.y.all)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Motif from all TFBS, to removing 1kb from promoters, to removing 10kb from promoters")
print(m.zscore.arrows.10kb)

# what's the average change? 

pdf("~/Dropbox/scCHiC_figs/FIG4_BM/analyses/zscore_promoter_enhancer_analysis.pdf", useDingbats = FALSE)
print(m.zscore.arrows)
print(m.zscore.arrows.10kb)
dev.off()

# Cluster TFs by distances in cell space ----------------------------------

# merge the matrices then do UMAP on the terms

act.mat.merged <- spread(act.long.merged, key = "cell", value = "activity")

# UMAP on the terms
motifs.filt2 <- subset(zscores.merged, !is.na(motif.lab))$motif.lab
motifs.filt <- zscores.merged$motif
jsub <- subset(act.mat.merged, motif %in% motifs.filt)
# umap.out <- umap(t(subset(jsub, select = -motif)))  # maybe interesting? need to color by batch
umap.out <- umap(subset(jsub, select = -motif))

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
umap.out.long <- as.data.frame(umap.out$layout)
umap.out.long$motif <- jsub$motif
umap.out.long$motif.lab <- sapply(umap.out.long$motif, function(x) ifelse(x %in% motifs.filt2, x, ""))
umap.out.long <- left_join(umap.out.long, subset(zscores.merged.mean, select = c(motif, zscore.mean)))

m.tfclstr <- ggplot(umap.out.long, aes(x = V1, y = V2, label = motif.lab, colour = zscore.mean)) + 
  geom_point() + geom_text_repel() + 
  xlab("Umap1 Cell-Cell Dist") + ylab("Umap2 Cell-Cell Dist") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient()
  # scale_color_gradient2(low = muted("blue"), mid = "grey80", high = muted("red"), midpoint = mean(umap.out.long$zscore.mean)) 



# Plot examples  ----------------------------------------------------------

nn.lst <- c(40, 35)
jmetric='euclidean' 
jmindist=0.4
jseed=123
custom.settings.lst <- lapply(nn.lst, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))

# Load LDA Output
dirmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"

inf.lst <- lapply(jchips, function(jchip){
  inf <- file.path(dirmain, paste0("lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_10_15_20_25.Robj"))
  return(inf)
})

out.lda.lst <- lapply(inf.lst, LoadLDA)

# umap.lda.lst <- lapply(out.lda.lst, DoUmapFromLDA, custom.settings)
umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

names(umap.lda.lst) <- jchips
names(mara.outs) <- jchips

dat.merged.lst <- lapply(jchips, function(jchip) left_join(umap.lda.lst[[jchip]], mara.outs[[jchip]]$act.long))
names(dat.merged.lst) <- jchips

# join umaplong with actlong
# jmotif <- "Cebpb"
# jmotif <- "Spib"

pdf(plotf, useDingbats = FALSE)
print(m.zscore)
print(m.tfclstr)
# plot the UMAPs 
for (jmotif in motifs.filt2){
  plts.lst <- lapply(jchips, function(jchip){
    return(PlotMotifInUmap(jmotif, dat.merged.lst[[jchip]], mara.outs[[jchip]]$zscores))
  })
  multiplot(plts.lst[[1]], plts.lst[[2]], cols = 2)
}
dev.off()