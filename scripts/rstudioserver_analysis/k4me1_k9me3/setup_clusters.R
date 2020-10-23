# Jake Yeung
# Date of Creation: 2020-10-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/k4me1_k9me3/setup_clusters.R
# set up clusters


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDAs  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl")
dir.create(outdir)
jmark1 <- "H3K4me1"
jmark2 <- "H3K9me3"
jmarkdbl <- "H3K4me1xH3K9me3"
jmarks <- c(jmark1, jmark2); names(jmarks) <- jmarks
outname1 <- paste0("cluster_tables_", jmark1, "_BM_all_round2.txt")
outname2 <- paste0("cluster_tables_", jmark2, "_BM_all_round2.txt")
outnamedbl <- paste0("cluster_tables_", jmarkdbl, "_BM_all_round2.txt")
outf1 <- file.path(outdir, outname1)
outf2 <- file.path(outdir, outname2)
outfdbl <- file.path(outdir, outnamedbl) 
outpdf <- paste0("plots_", jmark1, "_", jmark2, "_BM_all_round2.pdf")

pdf(outpdf, useDingbats = FALSE)



# inf1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_varfilt/lda_outputs.count_mat.", jmark1, ".varfilt_1.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1.varfilt_1.5.K-30.Robj"))
# assertthat::assert_that(file.exists(inf1))
# inf2 <- file.path(hubprefix, paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_varfilt/lda_outputs.count_mat.", jmark2, ".varfilt_1.5.K-30.binarize.FALSE/ldaOut.count_mat.H3K9me3.varfilt_1.5.K-30.Robj"))
# assertthat::assert_that(file.exists(inf2))

lda.outs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_varfilt/lda_outputs.count_mat.", jmark, ".varfilt_1.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".varfilt_1.5.K-30.Robj"))
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result, jsep = "")
  
  # do LDA
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  return(dat.umap)
})

tm.lsts <- lapply(jmarks, function(jmark){
    print(jmark)
    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_varfilt/lda_outputs.count_mat.", jmark, ".varfilt_1.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".varfilt_1.5.K-30.Robj"))
    assertthat::assert_that(file.exists(inf))
    load(inf, v=T)
    
    tm.result <- posterior(out.lda)
    tm.result <- AddTopicToTmResult(tm.result, jsep = "")
    return(tm.result)
})


# Check louvains ----------------------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(lda.outs[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() +
    scale_color_manual(values = cbPalette) + 
    theme_bw() +
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
})


# Make pretty H3K4me1 -----------------------------------------------------


neutros <- c("3", "12", "2")
dends <- c("9", "10")
eryths <- c("8", "7")
bcells1 <- c("6", "5")
bcells2 <- c("14")
nk <- c("11")
pdends <- c("13")
basos <- c("1")
hspcs <- c("4")

# remove mystery 
badlouvs <- c("15")

dat.umap1.clean <- lda.outs[[jmark1]] %>%
  # filter(!louvain %in% bad.louvs) %>%
  mutate(louv.orig = louvain,
         louvain = ifelse(louvain %in% neutros, "Neutrophils", louvain),
         louvain = ifelse(louvain %in% dends, "DCs", louvain),
         louvain = ifelse(louvain %in% eryths, "Erythroblasts", louvain),
         louvain = ifelse(louvain %in% nk, "NKs", louvain),
         louvain = ifelse(louvain %in% pdends, "pDCs", louvain),
         louvain = ifelse(louvain %in% basos, "Basophils", louvain),
         louvain = ifelse(louvain %in% hspcs, "HSPCs", louvain),
         louvain = ifelse(louvain %in% bcells2, "Bcells2", louvain),
         louvain = ifelse(louvain %in% bcells1, "Bcells1", louvain),
         louvain = ifelse(louvain %in% badlouvs, NA, louvain),
         cluster = louvain,
         louvain = louv.orig)

ggplot(lda.outs[[jmark1]], aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  ggtitle(jmark1) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


m.clean1 <- ggplot(dat.umap1.clean, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  ggtitle(jmark1) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.clean1)



# Check H3K9me3 -----------------------------------------------------------

ggplot(lda.outs[[jmark2]], aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# just call it louvainX. Remove bad louvains

badlouvs <- c("louvain3", "louvain4")

dat.umap2.clean <- lda.outs[[jmark2]] %>%
  rowwise() %>%
  mutate(louv.orig = louvain,
         louvain = paste("louvain", louvain, sep = ""), 
         louvain = ifelse(louvain %in% badlouvs, NA, louvain),
         cluster = louvain)

m.clean2 <- ggplot(dat.umap2.clean, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  ggtitle(jmark2) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.clean2)


# Check dbl ---------------------------------------------------------------

inf.dbl <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K4me1xH3K9me3.match_dbl.K-30.Robj"
load(inf.dbl, v=T)

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result, jsep = "")
# do LDA
dat.umap.dbl <- DoUmapAndLouvain(tm.result$topics, jsettings)
lda.outs[[jmarkdbl]] <- dat.umap.dbl

ggplot(lda.outs[[jmarkdbl]], aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umapdbl.clean <- lda.outs[[jmarkdbl]] %>%
  mutate(louv.orig = louvain,
         cluster = "louvainX")

m.cleandbl <- ggplot(dat.umapdbl.clean, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  ggtitle(jmarkdbl) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.cleandbl)


dev.off()

# 
# # Check variane -----------------------------------------------------------
# 
# dat.impute.log <- log2(t(tm.lsts$H3K9me3$topics %*% tm.lsts$H3K9me3$terms))
# jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# dat.var <- CalculateVarAll(dat.impute.log, jchromos)
# 
# 
# dat.merge <- left_join(dat.umap2.clean, dat.var)
# 
# ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
#   geom_point() + 
#   scale_color_viridis_c(direction = -1) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Write outputs -----------------------------------------------------------



fwrite(x = dat.umap1.clean %>% dplyr::select(c(cell, cluster, umap1, umap2, louvain)), file = outf1, sep = "\t", na = NA, quote = FALSE)
fwrite(x = dat.umap2.clean %>% dplyr::select(c(cell, cluster, umap1, umap2, louvain)), file = outf2, sep = "\t", na = NA, quote = FALSE)
fwrite(x = dat.umapdbl.clean %>% dplyr::select(c(cell, cluster, umap1, umap2, louvain)), file = outfdbl, sep = "\t", na = NA, quote = FALSE)

