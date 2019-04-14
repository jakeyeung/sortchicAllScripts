# Jake Yeung
# Date of Creation: 2019-03-27
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_TFactivity_LDAoutput_Annotations_Rdata_build95_allmarks.R
# All marks


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

library(igraph)

library(GGally)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaData.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")


# Functions ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3")
jmarks.repress <- c("H3K27me3", "H3K9me3")
jmarks.all <- c(jmarks, jmarks.repress)
jmarks.sub <- jmarks

# jmarks.all <- jmarks.sub
names(jmarks.all) <- jmarks.all


barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))



# Load LDAs ---------------------------------------------------------------
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
# out.lda.lst <- out.lda.new.lst[c("H3K4me1", "H3K4me3")]
out.lda.lst <- out.lda.new.lst

infs <- list(inf1, inf2, inf3, inf4)
names(infs) <- c(jmarks, jmarks.repress)
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark.all, inf=inf, convert.chr20.21.to.X.Y = TRUE), jmarks.all, infs, SIMPLIFY = FALSE)

tm.result.lst <- lapply(out.lda.new.lst, function(x) posterior(x))
count.imputed.lst <- lapply(tm.result.lst, function(tm.result) log10(t(tm.result$topic %*% tm.result$terms)))

annots.lst <- lapply(c(jmarks, jmarks.all), function(jmark) out.objs[[jmark]]$regions.annot)


# Get dirs ----------------------------------------------------------------

# match mark-mouse-repS{1-9} -> mark-mouse-rep{1,2}
# H3K4me1.BM.m1.S9.AH3VGVBGX9.TCACATGT -> H3K4me1.BM.m1.rep1 ... 


jsuffix <- "build95"
jdist <- "0"


marabase <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_build95.cells_from_bin_analysis"
msuffix <- "K50"
# msuffix <- "Kbest"
# msuffix <- "Kbest_techrepnamefix"
# msuffix <- "Kbest_techrepnamefix"

mdirs <- lapply(jmarks.sub, function(jmark){
  mdir <- file.path(marabase,
                    paste0("hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE_", msuffix, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--", msuffix),
                    paste0("hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE_", msuffix))
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- file.path(marabase,
                    paste0("hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE_", msuffix, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--", msuffix),
                    paste0("hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE_", msuffix))
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})



mdirs <- c(mdirs, mdirs.repress)
# mdirs <- c(mdirs)


switch.rep.hash <- GetRepSwitchHash(experihash)
mara.outs <- lapply(mdirs, LoadMARA, fix.tech.rep = FALSE, rep.prefix = "", swap.tech.rep = switch.rep.hash)
names(mara.outs) <- jmarks.all

# get colnames


# test <- LoadMARA(mdirs[[1]], fix.tech.rep = FALSE, rep.prefix = "rep", swap.tech.rep = TRUE)

# fix replicate name S9 -> rep2 for example

names(mara.outs) <- jmarks.all

head(mara.outs[[1]]$act.long)

act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long, mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)


# act.long.merged <- rbind(mara.outs[[1]]$act.long)
zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, mara.outs[[3]]$zscores, mara.outs[[4]]$zscores), left_join, by = "motif")
# zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores), left_join, by = "motif")
cnames <- c("motif", paste("zscore", jmarks.all, sep = "."))

colnames(zscores.merged) <- cnames

zscore.thres <- 0.75
zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
  ifelse(max(jrow[[2]], jrow[[3]], jrow[[4]], jrow[[5]]) > zscore.thres, jrow[[1]], NA)
  # ifelse(max(jrow[[2]]) > zscore.thres, jrow[[1]], NA)
})

zscores.merged.mean <- zscores.merged
zscores.merged.mean$zscore.mean <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, mean)
zscores.merged.mean$zscore.max <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, max)

zscores.merged.mean <- zscores.merged.mean %>%
  arrange(desc(zscore.mean))

# Show top hits -----------------------------------------------------------




# Make umaps --------------------------------------------------------------


jmetric.louv='euclidean' 
jmindist.louv=0.3
jseed.louv=123

nn.louv.new <- c(28, 35, 33, 31)
# nn.louv.new <- c(51, 51, 33, 31)
jmindist.new <- c(0.2, 0.15, 0.3, 0.3)
# nn.new <- c(40, 30, 45, 27)
nn.new <- c(48, 15, 45, 27)
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
  # prefactor <- ifelse(jmark %in% c("H3K4me1", "H3K4me3"), -1, 1)  # flip axis
  dat.umap <- dat.umap.new.lst[[jmark]]
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr.hash.new.lst[[jmark]][[x]]))
  print(head(dat.umap.long))
  dat.umap.long$mark <- jmark
  return(dat.umap.long)
})


# do UMAP then merge with activities
custom.settings.lst <- custom.settings.new.lst
umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)
names(umap.lda.lst) <- jmarks.all

# jmark <- "H3K4me3"
dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

print(head(dat.merged.lst))


# Plot --------------------------------------------------------------------

jmotif <- "Mecp2"
jmotif <- "Ebf3"
jmotif <- "Cebpd"
jmotif <- "Gata3"
jmotif <- "Brca1"
jmotif <- "Trp63"
jmotif <- "Brca1"
jmotif <- "Tfdp1"

jmotif <- "Tal1"
jmotif <- "Sox6"
jmotif <- "Brca1"
jmotif <- "Brca1"
jmotif <- "Tfdp1"
jmotif <- "Tal1"
jmotif <- "Hoxc6"
jmotif <- "Six1"
jmotif <- "Stat6"
jmotif <- "Zfp281"
jmotif <- "Ezh2"

jmotif <- "Ebf3"
jmotif <- "Snai1"
jmotif <- "Bptf"
jmotif <- "Pml"
jmotif <- "Ebf1"
jmotif <- "Cebpb"

jmotif <- "Hoxc6"

jmotif <- "Stat6"

jcolvec <- c("grey80", "grey50", "darkblue")
plts.lst <- lapply(jmarks.all, function(jmark){
  return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 2, colvec = jcolvec))
})

multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)

# multiplot(plts.lst[[1]], cols = 2)


# Save objects ------------------------------------------------------------

# # save data in best way
# jmark.keep <- c("H3K4me1")
# 
# dat.merged.lst.H3K4me1_only <- dat.merged.lst[jmark.keep]
# mara.outs.H3K4me1_only <- mara.outs[jmark.keep]
system.time(
  save(dat.umap.long.new.lst, dat.merged.lst, mara.outs, custom.settings.new.lst, tm.result.lst, count.imputed.lst, annots.lst, count.mat.lst, 
       file = paste0("~/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.", Sys.Date(), ".RData"))
)


