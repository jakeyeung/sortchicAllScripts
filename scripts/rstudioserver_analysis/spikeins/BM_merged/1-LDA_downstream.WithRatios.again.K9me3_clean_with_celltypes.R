# Jake Yeung
# Date of Creation: 2020-10-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.celltyping.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}



# Functions ---------------------------------------------------------------




# Constnats ---------------------------------------------------------------


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all"
dir.create(outdir)

# inmain <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046")
inmain <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all")
# load(inf.lda, v=T)

jsuffix <- ".filt_0.15_0.95_counts_and_l2r.K-30"


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

outs.lst <- lapply(jmarks, function(jmark){
  # dname <- paste0("lda_outputs.count_mat_", jmark, jsuffix, ".binarize.FALSE")
  dname <- paste0("lda_outputs.count_mat.", jmark, jsuffix, ".binarize.FALSE")
  indir <- file.path(inmain, dname)
  # fname <- paste0("ldaOut.count_mat_", jmark, "_counts_filt.2020-09-12.K-30.Robj")
  fname <- paste0("ldaOut.count_mat.", jmark, jsuffix, ".Robj")
  inf.lda <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.lda))
  print(inf.lda)
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  dat.umap$mark <- jmark
  return(list(dat.umap = dat.umap, tm.result = tm.result))
})

dat.umaps.lst <- lapply(outs.lst, function(jlst){
  return(jlst$dat.umap)
})

tm.result.lst <- lapply(outs.lst, function(jlst){
  return(jlst$tm.result)
})



# Load K9me3 separately cleaned up  ---------------------------------------

jmarkrep <- "H3K9me3"
inf.k9me3 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_H3K9me3_badclustremoved/lda_outputs.countmat_H3K9me3_newonly_badclustremoved.K-30.binarize.FALSE/ldaOut.countmat_H3K9me3_newonly_badclustremoved.K-30.Robj"
load(inf.k9me3, v=T)

tm.result.k9me3 <- posterior(out.lda)
dat.umap.k9me3 <- DoUmapAndLouvain(tm.result.k9me3$topics, jsettings)
dat.umap.k9me3$mark <- jmarkrep

ggplot(dat.umap.k9me3, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tm.result.lst[[jmarkrep]] <- tm.result.k9me3
dat.umaps.lst[[jmarkrep]] <- dat.umap.k9me3


jmarksall <- c(jmarks, jmarkrep)

for (jmark in jmarksall){
  m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}


dat.umaps.long <- bind_rows(dat.umaps.lst) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         jrep = GetRepBM(experiname = experi), 
         stype = AnnotateSortFromLayoutBMall(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jrep = jrep, jmark = mark))

dat.meta <- subset(dat.umaps.long, select = c(cell, experi, plate, rowcoord, colcoord, stype))

# Variance stuff ----------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var.lst <- lapply(jmarks, function(jmark){
  tm.result <- tm.result.lst[[jmark]]
  dat.impute <- t(log2(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute, jchromos)
}) 

dat.var <- dat.var.lst %>%
  bind_rows()

dat.merge <- left_join(dat.umaps.long, dat.var)


# Add log2ratio -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/spikein_info_BM_round2_all.txt"

dat.spikeins.mat <- fread(inf.spikeins) %>%
  rowwise() %>%
  mutate(plate = as.numeric(SplitGetLast(experi, jsplit = "-")),
         jrep = GetRepBM(experiname = experi), 
         stype = AnnotateSortFromLayoutBMall(plate = plate, rowcoord = rowcoord, colcoord = colcoord, jrep = jrep, jmark = mark))


jsub <- subset(dat.spikeins.mat, select = c(samp, chromocounts, spikeincounts)) %>%
  rowwise() %>%
  mutate(l2r = log2(chromocounts / spikeincounts))

dat.merge <- left_join(dat.merge, jsub, by = c("cell" = "samp"))


# Plot umap with log2ratio  -----------------------------------------------

ggplot(dat.merge, aes(x = umap1, y = umap2, color = l2r)) + 
  geom_point() + 
  facet_wrap(~mark) + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show boxplots
ggplot(dat.merge, aes(x = stype, y = l2r)) + 
  geom_boxplot() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# do K9me3 only 

louv2clst <- list("1" = "HSPCs", 
                  "2" = "Myeloid", 
                  "3" = "Myeloid",
                  "4" = "Eryth",
                  "5" = "Lymphoid")

jhash <- hash::hash(louv2clst)

dat.k9me3 <- subset(dat.merge, mark == "H3K9me3") %>%
  rowwise() %>%
  mutate(cluster = AssignHash(x = as.character(louvain), louv2clst, null.fill = NA))

ggplot(dat.k9me3, aes(x = cluster, y = l2r)) + 
  geom_boxplot() + 
  theme_bw() + 
  facet_wrap(~mark) + 
  ylab("log2(chromocounts/spikeins)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
