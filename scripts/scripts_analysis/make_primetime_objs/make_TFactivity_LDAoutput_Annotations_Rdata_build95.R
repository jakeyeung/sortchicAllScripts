# Jake Yeung
# Date of Creation: 2019-03-26
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_TFactivity_LDAoutput_Annotations_Rdata_build95.R
# description# 

# rm(list=ls())

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
source("scripts/Rfunctions/GetMetaData.R")


# Functions ---------------------------------------------------------------




# Get dirs ----------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3")
jsuffix <- "build95"
jdist <- "0"

marabase <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_build95"
msuffix <- "K50"

mdirs <- lapply(jmarks, function(jmark){
  mdir <- file.path(marabase,
                    paste0("hiddenDomains_cellmin_550-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE_", msuffix, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--", msuffix),
                    paste0("hiddenDomains_cellmin_550-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE_", msuffix))
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.repress <- c("H3K27me3", "H3K9me3")
# mdirs.repress <- lapply(jmarks.repress, function(jmark){
#   mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
#                  jmark, 
#                  ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
#                  "/",
#                  "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
#   print(mdir)
#   assertthat::assert_that(dir.exists(mdir))
#   return(mdir)
# })

# jmarks.all <- c(jmarks, jmarks.repress)
jmarks.all <- c(jmarks)
names(jmarks.all) <- jmarks.all

# mdirs <- c(mdirs, mdirs.repress)
mdirs <- c(mdirs)
mara.outs <- lapply(mdirs, LoadMARA, fix.tech.rep = TRUE)

# fix replicate name S9 -> rep2 for example

names(mara.outs) <- jmarks.all

head(mara.outs[[1]]$act.long)

# act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long, mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)
act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long)
# zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, mara.outs[[3]]$zscores, mara.outs[[4]]$zscores), left_join, by = "motif")
zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores), left_join, by = "motif")
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

print(m.zscore)

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

# Plot first two motifs label to hits

m.two <- ggplot(zscores.merged, aes(x = zscore.H3K4me1, y = zscore.H3K4me3, label = motif.lab)) + 
  geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlim(c(0, 3)) + ylim(c(0, 3))
print(m.two)



# Load build95 LDA --------------------------------------------------------


# jmark <- "H3K4me1"
out.lda.new.lst <- list()
inf1 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-15_20_25_30_35.Robj")
load(inf1, v=T)
out.lda.new.lst[["H3K4me1"]] <- ChooseBestLDA(out.lda)
inf2 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-15_20_25_30_35.Robj")
load(inf2, v=T)
out.lda.new.lst[["H3K4me3"]] <- ChooseBestLDA(out.lda)
inf3 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj")
load(inf3, v=T)
out.lda.new.lst[["H3K27me3"]] <- ChooseBestLDA(out.lda)
inf4 <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
load(inf4, v=T)
out.lda.new.lst[["H3K9me3"]] <- ChooseBestLDA(out.lda)

topics.mat.new.lst <- lapply(out.lda.new.lst, function(out.lda){
  return(posterior(out.lda)$topics)
})

infs.nobin <- list("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
names(infs.nobin) <- c(jmarks, jmarks.repress)

out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs.nobin, SIMPLIFY = FALSE)
count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

# waiting for repressive calculations to finish
out.lda.lst <- out.lda.new.lst[c("H3K4me1", "H3K4me3")]

infs <- list(inf1, inf2, inf3, inf4)
names(infs) <- c(jmarks, jmarks.repress)
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf, convert.chr20.21.to.X.Y = TRUE), jmarks.all, infs, SIMPLIFY = FALSE)

tm.result.lst <- lapply(out.lda.new.lst, function(x) posterior(x))
count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topic %*% tm.result$terms)))

annots.lst <- lapply(c(jmarks, jmarks.all), function(jmark) out.objs[[jmark]]$regions.annot)


# Make umaps --------------------------------------------------------------


jmetric.louv='euclidean' 
jmindist.louv=0.3
jseed.louv=123

nn.louv.new <- c(28, 35, 33, 31)
jmindist.new <- c(0.2, 0.15, 0.3, 0.3)
nn.new <- c(40, 30, 45, 27)
# custom.settings.new.lst <- lapply(nn.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))
custom.settings.new.lst <- mapply(function(x, y) GetUmapSettings(x, jmetric.louv, y, jseed.louv), nn.new, jmindist.new, SIMPLIFY = FALSE)
custom.settings.louv.new.lst <- lapply(nn.louv.new, function(x) GetUmapSettings(x, jmetric.louv, jmindist.louv, jseed.louv))

names(custom.settings.new.lst) <- c(jmarks, jmarks.repress)
names(custom.settings.louv.new.lst) <- c(jmarks, jmarks.repress)

# do UMAP on new settings
dat.umap.new.lst <- mapply(function(topics.mat, custom.settings) umap(topics.mat, config = custom.settings), 
                           topics.mat.new.lst, custom.settings.new.lst,
                           SIMPLIFY = FALSE)
# do Louvain on new.louv settings
clstr.hash.new.lst <- mapply(function(topics.mat, custom.settings) DoLouvain(topics.mat, custom.settings, dat.umap.long = NULL), 
                             topics.mat.new.lst, custom.settings.louv.new.lst)

#  assign cluster to umap
dat.umap.long.new.lst <- lapply(jmarks.all, function(jmark){
  dat.umap <- dat.umap.new.lst[[jmark]]
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr.hash.new.lst[[jmark]][[x]]))
  print(head(dat.umap.long))
  dat.umap.long$mark <- jmark
  return(dat.umap.long)
})

# do UMAP then merge with activities
custom.settings.lst <- custom.settings.new.lst[c("H3K4me1", "H3K4me3")]
umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

# Plot examples -----------------------------------------------------------

jcolvec <- c("gray80", "gray50", "darkblue")
jmotif <- "Cebpb"
jmotif <- "Pml"
jmotif <- "Ebf1"
jmotif <- "Gata3"
jmotif <- "Pml"

jmotif <- "Cebpa"

jmotif <- "Sox6"

jmotif <- "Ebf1"
jmotif <- "Tal1"
jmotif <- "Isl2"
jmotif <- "Atf2"
jmotif <- "Nfatc3"

jmotif <- "Pml"

plts.lst <- lapply(jmarks.all, function(jmark){
  return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec))
})
# multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)
multiplot(plts.lst[[1]], plts.lst[[2]], cols = 2)


# Save objects ------------------------------------------------------------

# save data in best way
jmark.keep <- c("H3K4me1")

dat.merged.lst.H3K4me1_only <- dat.merged.lst[jmark.keep]
mara.outs.H3K4me1_only <- mara.outs[jmark.keep]
system.time(
  save(dat.merged.lst.H3K4me1_only, mara.outs.H3K4me1_only, custom.settings.new.lst, tm.result.lst, count.imputed.lst, annots.lst, count.mat.lst, file = "~/data/scchic/robjs/TFactivity_genelevels_objects_build95.H3K4me1_activities_only.with_count_mats.RData")
)


# Save sample names  ------------------------------------------------------

# take bin analysis samples and convert to bam files



# 
# rm(list=ls())
# 
# 
# # Plot examples  ----------------------------------------------------------
# 
# nn.vec <- c(40, 35, 40, 40)
# jmindist.vec <- c(0.2, 0.1, 0.2, 0.1)
# jmetric <- "euclidean"
# jseed=123
# custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))
# 
# meanfilt <- 10
# 
# Kstr.bin <- "15_20_25_30_35"
# Kstr.nobin <- "15_20_25_30"
# 
# infs.nobin <- lapply(jmarks.all, function(jmark){
#   inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
#                 ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
#                 "/lda_out_meanfilt.BM-", jmark, 
#                 ".CountThres0.K-", Kstr.nobin, ".Robj")
#   assertthat::assert_that(file.exists(inf))
#   return(inf)
# })
# infs.bin <- lapply(jmarks.all, function(jmark){
#   inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
#                 ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
#                 "/lda_out_meanfilt.BM-", jmark, 
#                 ".CountThres0.K-", Kstr.bin, ".Robj")
#   assertthat::assert_that(file.exists(inf))
#   return(inf)
# })
# 
# infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
# # out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks.all, infs, SIMPLIFY = FALSE)
# 
# out.lda.lst <- lapply(infs, LoadLDA)
# 
# # out.lda.lst <- out.objs
# 
# umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)
# umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)
# 
# names(umap.lda.lst) <- jmarks.all
# names(mara.outs) <- jmarks.all
# 
# dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
# names(dat.merged.lst) <- jmarks.all
# 
# out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs, SIMPLIFY = FALSE)
# out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf), jmarks.all, infs.nobin, SIMPLIFY = FALSE)
# 
# # out.lda.lst <- out.objs
# 
# umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)
#  
# names(umap.lda.lst) <- jmarks.all
# names(mara.outs) <- jmarks.all
# 
# dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
# names(dat.merged.lst) <- jmarks.all
# 
# # tm.result.lst <- lapply(out.lda.lst, function(x) posterior(x))
# # count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topic %*% tm.result$terms)))
# count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))
# 
# # tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))
# 
# count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topics %*% tm.result$terms)))
# 
# annots.lst <- lapply(jmarks.all, function(jmark) out.objs[[jmark]]$regions.annot)
# 
# # # save data in best way 
# # system.time(
# #   save(dat.merged.lst, mara.outs, custom.settings.lst, count.mat.lst, tm.result.lst, count.imputed.lst, annots.lst, file = "~/data/scchic/robjs/TFactivity_genelevels_objects_build95.RData")
# # )
# 
# # rm(list=ls())
# 
# 
# 
# jcolvec <- c("gray80", "gray50", "darkblue")
# jmotif <- "Cebpb"
# plts.lst <- lapply(jmarks.all, function(jmark){
#   return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec))
# })
# multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)
# 
# 
# # plot a gene
# ref.mark <- "H3K4me1"
# jgene <- "S100a8"
# out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.jmark]])
# (jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]]))
# jscale.fac <- 1
# jpseudo <- 10^-6
# jsize <- 1
# PlotUmapAllMarks(jmarks.all, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, custom.settings.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)
# 
