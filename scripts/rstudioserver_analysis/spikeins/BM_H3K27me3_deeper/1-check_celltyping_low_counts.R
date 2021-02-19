# Jake Yeung
# Date of Creation: 2021-02-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_deeper/1-check_celltyping_low_counts.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# load metas --------------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(inf)
})

# Load LDA  ---------------------------------------------------------------


jmark <- "H3K27me3"

inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj")
load(inf, v=T)


tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(topics.mat = tm.result$topics, jsettings = jsettings)

dat.umap.annot <- left_join(dat.umap, subset(dat.metas[[jmark]], select = -c(umap1, umap2, louvain)))

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = log2(cuts_total))) + 
  geom_point() + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = log2(cuts_total))) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Get imputed -------------------------------------------------------------

dat.imputed <- t(tm.result$topics %*% tm.result$terms)

# avg across pseudbulk? 
cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
  
dat.imputed.sum <- SumAcrossClusters(dat.imputed, cnames.keep.lst)
dat.imputed.sum <- do.call(cbind, dat.imputed.sum)
dat.imputed.mean <- sweep(dat.imputed.sum, MARGIN = 2, STATS = colSums(dat.imputed.sum), FUN = "/")

probs.lst.filt <- as.list(as.data.frame(dat.imputed.mean))


# take raw
jsub <- (dat.metas[[jmark]] %>% arrange(cuts_total))
jcell <- jsub$cell[[1]]
count.vec <- count.mat[, jcell]

all.cells <- rownames(tm.result$topics); names(all.cells) <- all.cells

mfits <- FitMultinoms(count.filt = count.mat, all.cells = all.cells, probs.lst.filt = probs.lst.filt, exppower = 1)

LL.ctype.lst <- mfits

# Summarize ---------------------------------------------------------------

LL.dat <- SummarizeMultinomFits(LL.ctype.lst, count.mat, all.cells) %>%
  left_join(., dat.metas[[jmark]])


ggplot(LL.dat, aes(x = umap1, y = umap2, color = ctype.pred)) + 
  geom_point() +
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

LL.dat <- LL.dat %>%
  rowwise() %>%
  mutate(is.correct = ctype == ctype.pred)

head(LL.dat %>% arrange(cuts_total))

ggplot(LL.dat, aes(x = umap1, y = umap2, color = is.correct)) + 
  geom_point() +
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat, aes(x = log2(cuts_total), fill = is.correct)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

data.frame(print(head(LL.dat %>% arrange(cuts_total))))

LL.dat.sum <- LL.dat %>% 
  rowwise() %>%
  mutate(is.low = cuts_total < 1000) %>%
  group_by(is.low, is.correct) %>%
  summarise(ncorrect = length(is.correct)) %>%
  group_by(is.low) %>%
  mutate(fraccorrect = ncorrect / sum(ncorrect))

print(LL.dat.sum)

ggplot(LL.dat %>% mutate(is.low = cuts_total < 1000)) + 
  geom_point() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat.sum, aes(x = is.correct, y = fraccorrect, fill = is.low)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


save(LL.ctype.lst, LL.dat, dat.metas, count.mat, file = paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/multinomial_fits_celltype_single_cells/multinom_fits.", jmark, ".", Sys.Date(), ".RData"))

