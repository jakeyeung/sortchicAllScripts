# Jake Yeung
# Date of Creation: 2019-02-10
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/mara_downstream_both_marks.R
# Find motifs that are present in BOTH marks and which that are present in ONE mark

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

plotf <- "~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/mara/H3K4me1_and_H3K4me3_motifs_again.pdf"

jchips <- c("H3K4me1", "H3K4me3")

mdirs <- lapply(jchips, function(jchip){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jchip, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jchip, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

mara.outs <- lapply(mdirs, LoadMARA)

head(mara.outs[[1]]$act.long)


act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long)
zscores.merged <- left_join(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, by = "motif")
zscores.merged.mean <- zscores.merged %>%
  mutate(zscore.mean = (zscore.x + zscore.y) / 2)

zscore.thres <- 0.75
zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
  ifelse(max(jrow[[2]], jrow[[3]]) > zscore.thres, jrow[[1]], NA)
})

# ggplot(act.long.merged, aes(x = activity1, y = activity2, label = motif)) + 
#   geom_point() + 
#   geom_text()

m.zscore <- ggplot(zscores.merged, aes(x = zscore.x, y = zscore.y, label = motif.lab)) + 
  geom_point(alpha = 0.5) + theme_bw() + geom_text_repel() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab(paste("Motif Zscore:", jchips[[1]])) + ylab(paste("Motif Zscore:", jchips[[2]])) + 
  xlim(c(0, 2.5)) + ylim(c(0, 2.5))
print(m.zscore)

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

# plot(umap.out$layout[, 1], umap.out$layout[, 2], pch = 20)
# text(umap.out$layout[, 1], umap.out$layout[, 2], sapply(jsub$motif, function(x) ifelse(x %in% motifs.filt2, x, "")))


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