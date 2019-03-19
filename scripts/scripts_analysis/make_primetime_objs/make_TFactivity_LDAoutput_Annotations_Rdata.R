# Jake Yeung
# Date of Creation: 2019-03-19
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_TFactivity_LDAoutput_Annotations_Rdata.R
# Make TF activity and LDA output and annotations data for AvO
# 

rm(list=ls())

library(GGally)
library(purrr)

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



library(GGally)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")

# Functions ---------------------------------------------------------------




# Get dirs ----------------------------------------------------------------

plotf <- "~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/mara/all_4_marks_motifs_reordermotifs.pdf"

jmarks <- c("H3K4me1", "H3K4me3")

mdirs <- lapply(jmarks, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.repress <- c("H3K27me3", "H3K9me3")
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  print(mdir)
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all

mdirs <- c(mdirs, mdirs.repress)
mara.outs <- lapply(mdirs, LoadMARA)

head(mara.outs[[1]]$act.long)


act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long, mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)
# zscores.merged <- left_join(mara.outs[[3]]$zscores, mara.outs[[4]]$zscores, by = "motif")
# zscores.merged <- left_join(mara.outs[[3]]$zscores, mara.outs[[4]]$zscores, by = "motif")
zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, mara.outs[[3]]$zscores, mara.outs[[4]]$zscores), left_join, by = "motif")
cnames <- c("motif", paste("zscore", jmarks.all, sep = "."))

colnames(zscores.merged) <- cnames


zscore.thres <- 0.75
zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
  ifelse(max(jrow[[2]], jrow[[3]]) > zscore.thres, jrow[[1]], NA)
})

zscores.merged.mean <- zscores.merged
zscores.merged.mean$zscore.mean <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, mean)
zscores.merged.mean$zscore.max <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, max)

zscores.merged.mean <- zscores.merged.mean %>%
  arrange(desc(zscore.mean))

# pairs plot??
m.zscore <- ggpairs(zscores.merged, columns = paste("zscore", jmarks.all, sep = "."), 
                    lower = list(continuous = wrap("points", alpha = 0.2))) + theme_classic()



# Cluster TFs by distances in cell space ----------------------------------

# merge the matrices then do UMAP on the terms

act.mat.merged <- spread(act.long.merged, key = "cell", value = "activity")

# UMAP on the terms
motifs.filt2 <- subset(zscores.merged.mean, !is.na(motif.lab))$motif.lab
motifs.filt <- zscores.merged$motif
jsub <- subset(act.mat.merged, motif %in% motifs.filt)
# umap.out <- umap(t(subset(jsub, select = -motif)))  # maybe interesting? need to color by batch

# settings for UMAP
nn=10
jmetric='euclidean'
# jmetric='pearson2' 
jmindist=0.1
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
umap.out <- umap(subset(jsub, select = -motif), config = custom.settings)
# umap.out <- umap(jsub[, grepl("H3K27me3|H3K9me3", colnames(jsub))], config = custom.settings)
# umap.out <- umap(jsub[, grepl("H3K4me1|H3K27me3", colnames(jsub))], config = custom.settings)

# pca.out <- prcomp(subset(jsub, select = -motif), center = TRUE, scale. = TRUE)

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
umap.out.long <- as.data.frame(umap.out$layout)
# umap.out.long <- as.data.frame(pca.out$x[, 1:2]); colnames(umap.out.long) <- c("V1", "V2")
umap.out.long$motif <- jsub$motif
umap.out.long$motif.lab <- sapply(umap.out.long$motif, function(x) ifelse(x %in% motifs.filt2, x, ""))
umap.out.long <- left_join(umap.out.long, subset(zscores.merged.mean, select = c(motif, zscore.max, zscore.mean)))

m.tfclstr <- ggplot(umap.out.long, aes(x = V1, y = V2, label = motif.lab, colour = zscore.mean)) + 
  geom_point() + geom_text_repel() + 
  xlab("Umap1 Cell-Cell Dist") + ylab("Umap2 Cell-Cell Dist") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient()
print(m.tfclstr)


# Load LDA data -----------------------------------------------------------


# Plot examples  ----------------------------------------------------------

nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.1, 0.2, 0.1)
jmetric <- "euclidean"
jseed=123
custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
# out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks.all, infs, SIMPLIFY = FALSE)

out.lda.lst <- lapply(infs, LoadLDA)

# out.lda.lst <- out.objs

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

names(umap.lda.lst) <- jmarks.all
names(mara.outs) <- jmarks.all

dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs, SIMPLIFY = FALSE)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs.nobin, SIMPLIFY = FALSE)

# out.lda.lst <- out.objs

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

names(umap.lda.lst) <- jmarks.all
names(mara.outs) <- jmarks.all

dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

# tm.result.lst <- lapply(out.lda.lst, function(x) posterior(x))
# count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topic %*% tm.result$terms)))
count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

# tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))

count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topics %*% tm.result$terms)))

annots.lst <- lapply(jmarks.all, function(jmark) out.objs[[jmark]]$regions.annot)


# save data in best way 
system.time(
  save(dat.merged.lst, mara.outs, custom.settings.lst, count.mat.lst, tm.result.lst, count.imputed.lst, annots.lst, file = "~/data/scchic/robjs/TFactivity_genelevels_objects.RData")
)

# rm(list=ls())



jcolvec <- c("gray80", "gray50", "darkblue")
jmotif <- "Cebpb"
plts.lst <- lapply(jmarks.all, function(jmark){
  return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec))
})
multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)


# plot a gene
ref.mark <- "H3K4me1"
jgene <- "S100a8"
out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.jmark]])
(jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]]))
jscale.fac <- 1
jpseudo <- 10^-6
jsize <- 1
PlotUmapAllMarks(jmarks.all, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, custom.settings.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)

