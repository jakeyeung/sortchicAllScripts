# Jake Yeung
# Date of Creation: 2020-12-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/10-check_genes_associated_with_bins.R
# 

rm(list=ls())

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)
library(DescTools)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmarks <- c("H3K4me1", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks



keeptop <- 150
low.in.k9 <- TRUE
# low.in.k9 <- FALSE
# outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/heatmap_k9me3_k4me1_signif_bins_k9.highink9_", low.in.k9, ".", Sys.Date(), ".WithLogFCmaps.pdf")

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs"
intxt <- file.path(indir, paste0("correlate_k9_specific_bins.topn_", keeptop, ".LowInK9_", low.in.k9, ".2020-12-22.txt"))

intxt.k4 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/H3K9me3_analysis/heatmap_allmarks_signif_bins_k9.lowink9_TRUE.2020-12-23.AllInOne.params.H3K4me1.txt"
intxt.k9 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/H3K9me3_analysis/heatmap_allmarks_signif_bins_k9.lowink9_TRUE.2020-12-23.AllInOne.params.H3K9me3.txt"

dat.params.k4 <- fread(intxt.k4)
dat.params.k9 <- fread(intxt.k9)


dat.params <- fread(intxt)

bins.keep <- dat.params$bin
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.uniq <- unique(bins.keep)
dat.annot <- AnnotateCoordsFromList(coords.vec = bins.uniq, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

# annotate dat.params

subset(dat.annot$out2.df.closest, region_coord %in% jbins)$gene

dat.params.annot <- left_join(dat.params, subset(dat.annot$out2.df.closest, select = c(region_coord, gene)), by = c("bin" = "region_coord"))
dat.params.annot.k4 <- left_join(dat.params.k4, subset(dat.annot$out2.df.closest, select = c(region_coord, gene)), by = c("bin" = "region_coord"))
dat.params.annot.k9 <- left_join(dat.params.k9, subset(dat.annot$out2.df.closest, select = c(region_coord, gene)), by = c("bin" = "region_coord"))


jfilt <- subset(dat.params.annot, !is.na(gene))
jfilt.k4 <- left_join(subset(dat.params.annot.k4, !is.na(gene)), subset(jfilt, select = c(bin, label)))
jfilt.k9 <- left_join(subset(dat.params.annot.k9, !is.na(gene)), subset(jfilt, select = c(bin, label)))

# jfilt.merge <- left_join(jfilt.k4, jfilt.k9 by = c("bin"))

# get a eryth effect
jfilt.eryth <- subset(jfilt.k9, label == "Eryths") %>%
  arrange(Eryths.effect) %>%
  left_join(., jfilt.k4 %>% dplyr::rename(Eryths.effect.k4 = Eryths.effect) %>% dplyr::select(c(bin, Eryths.effect.k4)), by = c("bin"))
head(print(jfilt.eryth %>% dplyr::select(c(bin, Eryths.effect, Eryths.effect.k4, gene, label))))

library(ggrepel)
ggplot(jfilt.eryth, aes(x = Eryths.effect, y = Eryths.effect.k4, label = gene)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt.bcell <- subset(jfilt.k9, label == "Bcells") %>%
  arrange(Bcells.effect) %>%
  left_join(., jfilt.k4 %>% dplyr::rename(Bcells.effect.k4 = Bcells.effect) %>% dplyr::select(c(bin, Bcells.effect.k4)), by = c("bin")) 
head(print(jfilt.bcell %>% dplyr::select(c(bin, Bcells.effect, Bcells.effect.k4, gene, label))), n = 10)

ggplot(jfilt.bcell, aes(x = Bcells.effect, y = Bcells.effect.k4, label = gene)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt.granu <- subset(jfilt.k9, label == "Granulocytes") %>%
  arrange(Granulocytes.effect) %>%
  left_join(., jfilt.k4 %>% dplyr::rename(Granulocytes.effect.k4 = Granulocytes.effect) %>% dplyr::select(c(bin, Granulocytes.effect.k4)), by = c("bin"))
head(print(jfilt.granu %>% dplyr::select(c(bin, Granulocytes.effect, Granulocytes.effect.k4, gene, label))))

ggplot(jfilt.granu, aes(x = Granulocytes.effect, y = Granulocytes.effect.k4, label = gene)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt.hspcs <- subset(jfilt.k9, label == "HSPCs") %>%
  arrange(HSPCs.effect) %>%
  left_join(., jfilt.k4 %>% dplyr::rename(HSPCs.effect.k4 = HSPCs.effect) %>% dplyr::select(c(bin, HSPCs.effect.k4)), by = c("bin"))
head(print(jfilt.hspcs %>% dplyr::select(c(bin, HSPCs.effect, HSPCs.effect.k4, gene, label))))

ggplot(jfilt.hspcs, aes(x = HSPCs.effect, y = HSPCs.effect.k4, label = gene)) + 
  geom_point() + 
  geom_text_repel() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Find specific genes close to a gene  ------------------------------------



