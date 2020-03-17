# Jake Yeung
# Date of Creation: 2020-02-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/10-find_tx_diversity_in_stem_cells.R
# Order cells by "transcriptional diversity"? 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(ggrepel)



jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))


# Load GLMOUT -------------------------------------------------------------

inf <- "/home/jyeung/hpc/scChiC/from_rstudioserver/GLMPCA_outputs_copiedFromServer/PZ_H3K4me3.GLMPCA_var_correction.150.2020-02-01.binskeep_1000.devmode.RData"
load(inf, v=T)


tm.result <- list(topics = as.matrix(glm.out$factors), terms = t(as.matrix(glm.out$loadings)))
dat.impute.glm <- t(tm.result$topics %*% tm.result$terms)  #  already in log-odds

colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")


annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))


# Filter for genes only  --------------------------------------------------

jterms.keep.dat <- subset(annots.out$terms.annot, !is.na(gene))
jterms.keep <- subset(annots.out$terms.annot, !is.na(gene))$term

dat.keep <- dat.impute.glm[jterms.keep, ]

# calculate entropy for each gene

vec.entropy <- apply(dat.keep, 2, function(jcol) CalculateEntropy(exp(jcol), normalize.p = TRUE))
dat.entropy <- data.frame(cell = names(vec.entropy), S = vec.entropy, stringsAsFactors = FALSE)


# Do umap  ----------------------------------------------------------------


topics.mat <- tm.result$topics

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long) %>%
  left_join(., dat.entropy)

dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = plate)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~plate)

ggplot(dat.umap.long %>% mutate(S = ifelse(S < 18, NA, S)), aes(x = umap1, y = umap2, color = S)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = 1, na.value = "grey85")


# Fraction of reads in genes versus outside?  -----------------------------





