# Jake Yeung
# Date of Creation: 2020-12-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/19-make_gene_lists_for_filtering.R
# Make gene lists for filtering then rerun GLMPCA faster

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load gene lists  --------------------------------------------------------



# load H3K4me1 and find the Il6?

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.k4me1 <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.H3K4me1.dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.H3K4me1.dist_10000.K-30.Robj")
load(inf.k4me1, v=T)


tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load metadata find basophils  -------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/metadata_umap_celltype_cuts.H3K4me1.txt")
dat.meta <- fread(inf.meta)

ggplot(dat.meta, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.annot <- left_join(dat.umap, subset(dat.meta, select = c(cell, cluster)))

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~cluster) + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Order topics by entropy  ------------------------------------------------

topics.ordered <- OrderTopicsByEntropy(tm.result = tm.result)

jtopic <- "topic3"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1"
outpdf <- file.path(outdir, "H3K4me1_TSS_genesets.pdf")

pdf(outpdf, useDingbats = FALSE)

for (jtopic in topics.ordered$topic){
  jloadings <- data.frame(cell = rownames(tm.result$topics), loadings = tm.result$topics[, jtopic], stringsAsFactors = FALSE) %>%
    left_join(dat.meta, .)
  
  m <- ggplot(jloadings, aes(x = umap1, y = umap2, color = loadings)) + 
    geom_point() + 
    scale_color_viridis_c() + 
    theme_bw() + ggtitle(jtopic) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}
dev.off()

# get genes from topic20

jtopic <- "topic4"
top.genes <- sort(tm.result$terms[jtopic, ], decreasing = TRUE)
top.rnames <- names(rank(-sort(tm.result$terms[jtopic, ], decreasing = TRUE))[1:400])

jgenes <- sapply(top.rnames, function(x) strsplit(x, "\\.")[[1]][[4]])

# check giladi 

inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)

dat.sub <- subset(dat.public, gene %in% jgenes)

ggplot(dat.sub, aes(x = celltype, y = zscore)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Assign topics to genes  -------------------------------------------------

annot <- 
  list("topic30" = "pDCs",
       "topic3" = "NKs",
       "topic28" = "DCs", 
       "topic2" = "Bcells",
       "topic12" = "Bcells2",
       "topic16" = "Eryths", 
       "topic20" = "Basophils",
       "topic10" = "Granulocytes",
       "topic11" = "HSPCs",
       "topic4" = "HSPCs2")
  
# check
jnames <- names(annot)[order(unlist(annot))]
names(jnames) <- jnames

# get top genes and write

topn <- 400
dat.genes.k4me1 <- lapply(jnames, function(jname){
  clstr <- annot[[jname]]
  top.genes <- sort(tm.result$terms[jname, ], decreasing = TRUE)
  top.ranks <- rank(-sort(tm.result$terms[jname, ], decreasing = TRUE))[1:topn]
  top.rnames <- names(top.ranks)
  jgenes <- sapply(top.rnames, function(x) strsplit(x, "\\.")[[1]][[4]])
  dat.genes <- data.frame(jset = clstr, topic = jname, gene = top.rnames, symbol = jgenes, rnk = top.ranks, stringsAsFactors = FALSE)
}) %>%
  bind_rows()

# annotate colors 



outtxt <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.txt"
fwrite(dat.genes.k4me1, file = outtxt, sep = "\t")

inf.genes.old <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
dat.genes.old <- fread(inf.genes.old)
# filter out and add colors
jset2color <- hash::hash(unique(dat.genes.old$jset), unique(dat.genes.old$colorcode))
dat.genes.k4me1$colorcode <- sapply(dat.genes.k4me1$jset, function(x) AssignHash(x = x, jhash = jset2color, null.fill = NA))
dat.genes.sub <- subset(dat.genes.k4me1, !is.na(colorcode))
outf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt")
fwrite(x = dat.genes.sub, outf.genes)

# Check K4me1 bins  -------------------------------------------------------



