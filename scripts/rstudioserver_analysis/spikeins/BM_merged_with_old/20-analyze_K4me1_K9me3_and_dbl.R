# Jake Yeung
# Date of Creation: 2020-12-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/20-analyze_K4me1_K9me3_and_dbl.R
# Compare k4me1, k9me3 bins and also dbl (reseq)


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)


hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load raw cuts  ----------------------------------------------------------

indir.cuts <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables"))
infs.cuts <- list.files(indir.rz, pattern = ".*.binsize_50000.csv", full.names = TRUE)

indir.rz <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables"))
infs.rz <- list.files(indir.rz, pattern = ".*.RZ.csv", full.names = TRUE)

indir.chromo <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables"))
infs.chromo <- list.files(indir.chromo, pattern = ".*.NoChromo.csv", full.names = TRUE)

# Load cuts ---------------------------------------------------------------

dats.rz <- lapply(infs.rz, function(inf){
  ReadLH.SummarizeTA(inf)
}) %>%
  bind_rows()

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")
dats.chromo <- lapply(infs.chromo, function(infchromo){
  GetChromoCounts(infchromo) %>%
    filter(chromo == jspikeinchromo)
}) %>%
  bind_rows()



# Load bins ---------------------------------------------------------------

mat.lst <- lapply(infs.cuts, function(inf){
  ReadMatSlideWinFormat(inf)
}) 

rnames.all <- unique(unlist(lapply(mat.lst, rownames)))

mat.merge <- cbind.fill.lst(mat.lst, rnames.all)



# Check TA frac -----------------------------------------------------------

dat.meta.merge <- left_join(dats.rz, dats.chromo)

ggplot(dat.meta.merge, aes(x = log10(chromocounts), y = TA.frac)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.merge, aes(x = log2(chromocounts / spikeincounts), y = TA.frac)) + 
  geom_point(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Calculate fraction of bins with non-zero cuts ---------------------------

frac.nonzeros <- apply(mat.merge, MARGIN = 2, FUN = function(jcol){
  nnzero(jcol) / length(jcol)
})

dat.frac.nonzeros <- data.frame(samp = names(frac.nonzeros), frac.nzeros = frac.nonzeros, stringsAsFactors = FALSE)

dat.meta.merge2 <- left_join(dat.meta.merge, dat.frac.nonzeros)


ggplot(dat.meta.merge2, aes(x = log10(chromocounts), y = frac.nzeros)) +
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.merge2, aes(x = log2(chromocounts / spikeincounts), y = frac.nzeros)) +
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Filter cells ------------------------------------------------------------


l2r.min <- -1
ta.min <- 0.5
cuts.min <- 1000
nfrac.max <- 0.5


dat.meta.merge3 <- dat.meta.merge2 %>%
  rowwise() %>%
  mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
         colcoord = AddPlateCoordinates(samp)$colcoord,
         is.empty = rowcoord <= 8 & colcoord == 1,
         mark = "H3K4me1-H3K9me3",
         stype = ifelse(rowcoord <= 12, "Unenriched", "Linneg"),
         l2r = log2(chromocounts / spikeincounts), 
         is.good = !is.empty & l2r > l2r.min & TA.frac > ta.min & chromocounts > cuts.min & nfrac.max <= frac.nzeros)


# Keep good cells write to output for LDA  --------------------------------

cells.keep <- subset(dat.meta.merge3, !is.good)$samp

ckeep <- colnames(mat.merge) %in% cells.keep

mat.merge.filt <- mat.merge[, ckeep]

library(scchicFuncs)
library(irlba)
library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

lsi.out <- RunLSI(as.matrix(mat.merge.filt))

dat.umap.lsi <- DoUmapAndLouvain(lsi.out$u, jsettings)

ggplot(dat.umap.lsi, mapping = aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

outrds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.H3K4me1_H3K9me3.reseq/count_mat_cleaned_reseq.H3K4me1_H3K9me3.rds"

saveRDS(mat.merge.filt, file = outrds)




