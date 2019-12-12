# Jake Yeung
# Date of Creation: 2019-12-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/summarize_LDA_celltypes_all_marks.R
# All marks celltypes 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(ggrepel)

library(JFuncs)


knitr::opts_chunk$set(echo=F, warnings=F, message=T, fig.width=16, fig.asp=1, dpi=50, fig.align='center')

# Functions ---------------------------------------------------------------

PlotTopicLoadings <- function(jtop, jmark){
  # jtop <- "topic_11"
  col.i <- as.numeric(strsplit(jtop, "_")[[1]][[2]])
  
  jdat <- data.frame(topic.weight = dat.outputs[[jmark]]$tm.result$topics[, col.i], cell = rownames(dat.outputs[[jmark]]$tm.result$topics))
  
  dat.merge <- left_join(dat.outputs[[jmark]]$dat.merge, jdat)
  
  m.umap.topic <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = topic.weight)) + scale_color_viridis_c() + ggtitle(jmark, jtop) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m.umap.topic)
}


# Load LDA ----------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks




# Constants ---------------------------------------------------------------



jbin <- "FALSE"
jsuff <- "AllMerged"
# outdir <- "/home/jyeung/data/from_rstudioserver/scchic/pdfs/"
outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/LDA_downstream_ZF"

jwin <- 100000L

zscore.cutoff <- "top_500"

# load Chloe data
jchromos <- paste("chr", seq(25), sep = "")
inf.annot <- paste0("/home/jyeung/hpc/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed")
assertthat::assert_that(file.exists(inf.annot))


# normtype <- ""
normtype <- "_cpmnorm"
normtype2 <- "_nothrombo"
jdate0 <- ".2019-12-10"
jdate <- ".2019-12-10"

inf.WKM <- paste0("/home/jyeung/hpc/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, jdate0, ".rds")
assertthat::assert_that(file.exists(inf.WKM))

inf.WKMnothrombo <- paste0("/home/jyeung/hpc/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype2, jdate, ".rds")
assertthat::assert_that(file.exists(inf.WKMnothrombo))

# UMAP settings
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk$celltype <- as.factor(dat.bulk$celltype)

dat.bulk.nothrombo <- readRDS(inf.WKMnothrombo)
dat.bulk.nothrombo$celltype <- as.factor(dat.bulk.nothrombo$celltype)


outdir <- "~/data/from_rstudioserver/scchic/ZF_objs"
dir.create(outdir)
# load(file.path(outdir, "dat_across_marks_outputs.RData"), v=T)
load(file.path(outdir, "dat_across_marks_outputs.binactive_nobinrepress.RData"), v=T)

#' ## Show UMAPs across marks and conditions

m.umaps <- lapply(jmarks, function(jmark){
  dat.output <- dat.outputs[[jmark]]
  jsize <- 1
  m <- ggplot(dat.output$dat.merge %>% arrange(desc(cell.var.within.sum.norm)), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point(size = jsize) + scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + ggtitle(jmark)
  m.split <- ggplot(dat.output$dat.merge %>% arrange(desc(cell.var.within.sum.norm)), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point(size = jsize) + scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + facet_wrap(~experi) + 
    ggtitle(jmark)
  return(list(m = m, m.split = m.split))
})

#' In CD41low conditions, we see specific enrichment and depletions in different parts of the UMAP

lapply(list(m.umaps$H3K4me1$m, m.umaps$H3K4me3$m, m.umaps$H3K27me3$m, m.umaps$H3K9me3$m), function(x) return(x))

lapply(list(m.umaps$H3K4me1$m.split, m.umaps$H3K4me3$m.split, m.umaps$H3K9me3$m.split), function(x) return(x))
# JFuncs::multiplot(m.umaps$H3K4me1$m, m.umaps$H3K4me3$m, m.umaps$H3K27me3$m, m.umaps$H3K9me3$m, cols = 4)
# JFuncs::multiplot(m.umaps$H3K4me1$m.split, m.umaps$H3K4me3$m.split, m.umaps$H3K9me3$m.split, cols = 4)

#' ## Can we find HSC signature in UMAP? Let's do some analysis on H3K4me1

#' Here we see depletion of high variance cells at the top (erythyroblasts) and bottom left (monocytes/neutrophils), and an enrichment in cells at bottom (HSCs)

print(m.umaps$H3K4me1$m.split)

#' We find a topic corresponding to enrichment in Cd41-gated conditions (presumably HSCs)

keepn <- 100

# for bin=FALSE
# # jtop.hsc.h3k4me1 <- "topic_18"  # HSC only from Baron et al.
# jtop.hsc.h3k4me1 <- "topic_11"  # HSC based on both Baron et al. and Kobayashi et al.


# for bin=TRUE
jtop.hsc.h3k4me1 <- "topic_1"
m.umap.topic <- PlotTopicLoadings(jtop.hsc.h3k4me1, "H3K4me1")
print(m.umap.topic)

#' This topic is associated with HSC marker genes


jsub.terms <- subset(dat.outputs$H3K4me1$annot.pout$terms.annot, topic == jtop.hsc.h3k4me1) %>%
  dplyr::filter(rnk <= keepn)

m.top <- jsub.terms %>%
  ungroup() %>%
  mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
  ggplot(aes(x = term, y = log10(weight), label = gene)) +
  geom_point(size = 0.25) +
  theme_bw(8) +
  geom_text_repel(size = 2.5, segment.size = 0.1, segment.alpha = 0.25) +
  theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
  xlab("") + ylab("Log10 Bin Weight") +
  ggtitle(paste("Top peak weights for:", jtop.hsc.h3k4me1))

# check gene expression across genes
gfilt <- unique(jsub.terms$gene)
# gfilt <- bg.genes.all
# m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
  # mutate(celltype = ) %>%
  ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
  geom_point() + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle(jtop.hsc.h3k4me1)


#' Many genes in this topic have been previously reported as HSC markers

print(m.top)

#' runx1t1, meis1b, sox4a, gata2a (5 / 25) of these genes have been linked to recent ZF hematopoietic stem cell paper from Kobayashi et al. Scientific Reports 2019
#' 
#' https://www.nature.com/articles/s41598-019-50672-5


bg.genes.all <- unique(dat.outputs$H3K4me1$annot.pout$terms.annot$gene)
# gfilt <- sample(x = bg.genes.all, size = 25, replace = FALSE)


# Lymphoid markers (neg control)
# fg.published <- c("cd4-1", "cd8a", "gata3", "ighz", "igl4v8", "lck", "mhc2dab", "mhc2dbb", "pax5", "rag1", "traj15", "trdj2")
# HSC markers (what's interesting)

fg.published <- c("abcg2a", "apoeb", "cxcr4b", "egr1", "ENSDARG00000017320", "fgd5b", "fosb", "gata2a", "gata2b", "gfi1aa", "ikzf1", "kita", "meis1b", "mpl", "myb", "myca", "pbx1b", "runx1", "runx1t1", "sox4a", "tal1", "zfpm1")

fg.overlap <- intersect(gfilt, fg.published)
bg.overlap <- intersect(gfilt, bg.genes.all)

# Should be high enrichment
N.foverlap <- length(fg.overlap)
N.boverlap <- length(bg.overlap)

fish.mat <- matrix(data = c(N.foverlap, length(gfilt) - N.foverlap, N.boverlap, length(bg.genes.all) - N.boverlap), ncol = 2, nrow = 2, byrow = TRUE)
print(fish.mat)
jtest <- fisher.test(fish.mat)

print(jtest)  # 5/25 is highly significant

print(m.ctype)


#' Integrating with Baron et al. 2019 we see this set of genes is highly expressed in HSCs compared to other celltypes
#' (second place is thrombocytes, which seem to have some stem cell signature in my analysis)

#' Not sure how much of this is technical because in the pseudobulk analysis thrombocytes have an order of magnitude fewer total counts compared to other cell types. 




#' ## Can we reproduce this finding in H3K4me3?


print(m.umaps$H3K4me3$m.split)

# for bin=FALSE
# jtop.hsc.h3k4me3 <- "topic_29"  # HSC based on Baron et al. only
# jtop.hsc.h3k4me3 <- "topic_10"  # HSC topic based on both Baron et al. and Kobayashi et al.

# for bin=TRUE
jtop.hsc.h3k4me3 <- "topic_19"
jtop.hsc.h3k4me3 <- "topic_10"

#' UMAP is fuzzier in H3k4me3 versus H3K4me1 (H3K4me3 was less sequenced overall compared to H3K4me1)
m.umap.sc.h3k4me3 <- PlotTopicLoadings(jtop.hsc.h3k4me3, "H3K4me3")
print(m.umap.sc.h3k4me3)

jsub.terms <- subset(dat.outputs$H3K4me3$annot.pout$terms.annot, topic == jtop.hsc.h3k4me3) %>%
  dplyr::filter(rnk <= keepn)

m.top <- jsub.terms %>%
  ungroup() %>%
  mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
  ggplot(aes(x = term, y = log10(weight), label = gene)) +
  geom_point(size = 0.25) +
  theme_bw(8) +
  geom_text_repel(size = 2.5, segment.size = 0.1, segment.alpha = 0.25) +
  theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
  xlab("") + ylab("Log10 Bin Weight") +
  ggtitle(paste("Top peak weights for:", jtop.hsc.h3k4me3))

# check gene expression across genes
gfilt <- unique(jsub.terms$gene)
# m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
print(m.ctype)
# gfilt <- bg.genes.all
m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
  # mutate(celltype = ) %>%
  ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
  geom_point() + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle(jtop.hsc.h3k4me3)


fg.overlap <- intersect(gfilt, fg.published)
bg.overlap <- intersect(gfilt, bg.genes.all)

# Should be high enrichment
N.foverlap <- length(fg.overlap)
N.boverlap <- length(bg.overlap)

fish.mat <- matrix(data = c(N.foverlap, length(gfilt) - N.foverlap, N.boverlap, length(bg.genes.all) - N.boverlap), ncol = 2, nrow = 2, byrow = TRUE)
jtest <- fisher.test(fish.mat)

print(jtest)  # not significant? Possibly because H3K4me3 is much lower sequenced than H3K4me1. The stem cell markers genes are lower in the ranking of the topic weights

#' In Baron et al. this topic seems robustly highly expressed in HSPCs. Again, second place is thrombocytes and I'm not sure whether this is technical or not (10 times fewer reads in thrombocytes vs other celltypes)   

print(m.ctype)


#' ## Stem cells in H3K9me3

#' We don't have good markers for H3K9me3 for HSCs, but from our previous analyses, enrichment from the cd41-gated cells seem to have signatures for HSCs. Depletions are strongest for
#' erythrocytes and monocytes. So let's see whether we see this enrichment and depletion in the H3K9me3 UMAP

print(m.umaps$H3K9me3$m.split)

jtop.hsc.h3k9me3 <- "topic_15"

#' UMAP is fuzzier in H3k4me3 versus H3K4me1 (H3K4me3 was less sequenced overall compared to H3K4me1)
m.umap.sc.h3k9me3 <- PlotTopicLoadings(jtop.hsc.h3k9me3, "H3K9me3")
print(m.umap.sc.h3k9me3)


jsub.terms <- subset(dat.outputs$H3K9me3$annot.pout$terms.annot, topic == jtop.hsc.h3k9me3) %>%
  dplyr::filter(rnk <= keepn)

m.top <- jsub.terms %>%
  ungroup() %>%
  mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
  ggplot(aes(x = term, y = log10(weight), label = gene)) +
  geom_point(size = 0.25) +
  theme_bw(8) +
  geom_text_repel(size = 2.5, segment.size = 0.1, segment.alpha = 0.25) +
  theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
  xlab("") + ylab("Log10 Bin Weight") +
  ggtitle(paste("Top peak weights for:", jtop.hsc.h3k9me3))

# check gene expression across genes
gfilt <- unique(jsub.terms$gene)
# m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
# gfilt <- bg.genes.all
m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
  # mutate(celltype = ) %>%
  ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
  geom_point() + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle(jtop.hsc.h3k9me3)

## #+ fig.width=20, fig.asp=1, dpi=50, fig.align='center'

print(m.ctype)


#' ## K27me3 landscape in kidney marrow

#' K27me3 for cd41-gated cells data still not in yet. To be continued here...

print(m.umaps$H3K27me3$m.split)

#' ## Discussion




