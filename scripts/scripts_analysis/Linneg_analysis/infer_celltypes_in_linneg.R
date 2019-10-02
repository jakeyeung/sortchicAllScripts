# Jake Yeung
# Date of Creation: 2019-09-15
# File: ~/projects/scchic/scripts/scripts_analysis/Linneg_analysis/infer_celltypes_in_linneg.R
# Infer celltypes from linneg

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(umap)
library(here)

library(Matrix)

library(scchicFuncs)

library(ggrepel)

setwd(here())

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Functions ---------------------------------------------------------------



DoFishersTestCtype <- function(jctype, celltype.fc){
  conting.tab <- spread(celltype.fc %>% 
                          mutate(ctype.pred = ifelse(ctype.pred == jctype, jctype, "zNot")) %>% 
                          dplyr::select(batch, ctype.pred, cell.count) %>%
                          group_by(ctype.pred, batch) %>%
                          summarise(cell.count = sum(cell.count)), 
                        key = batch, value = cell.count) %>%
    as.data.frame() %>%
    ungroup()
  if (any(is.na(conting.tab))){
    warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
    return(NULL)
  }
  # rownames
  rownames(conting.tab) <- conting.tab$ctype.pred
  conting.tab$ctype.pred <- NULL
  
  # make odds ratio interpretable:( zNot_treat / ctype_treat ) / ( zNot_ctrl / ctype_ctrl )
  conting.tab <- t(conting.tab)
  # swap rows
  conting.tab <- conting.tab[c(2, 1), ]
  print(conting.tab)
  hyp.test <- fisher.test(x = conting.tab)
  return(hyp.test)
}


# Settings ----------------------------------------------------------------


jmark <- "H3K4me3"

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

custom.settings <- jsettings


# Constants for umap ------------------------------------------------------


mindist <- 0.4
nneigh <- 34
randomstate <- 123
jmetric <- "euclidean"

# umap.settings <- GetUmapSettings(nn = nneigh, jmindist = mindist, jmetric = jmetric, seed = randomstate)
umap.settings <- umap.defaults
umap.settings$min_dist <- mindist; umap.settings$n_neighbors <- nneigh

# Load original data  -----------------------------------------------------

# inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg_from_projection/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50_mindist_0.4_mindist_processed_lda.Rdata"
# inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.H3K4me3.RData"
# inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

# umap.orig <- umap(out.objs$tm.result$topics, config = umap.settings)
umap.orig <- umap(out.objs$tm.result$topics, config = custom.settings)

umap.orig.long <- data.table(umap1 = umap.orig$layout[, 1], umap2 = -1 * umap.orig$layout[, 2], samp = rownames(umap.orig$layout), batch = "Ctrl")

ggplot(umap.orig.long, aes(x = umap1, y = umap2)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# load new data
# inf.new <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg_from_projection/bin_PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.2019-06-17.RData"
inf.new <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg_from_projection/bin_TRUE_PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.from_louvain.RData"
load(inf.new, v=T)


umap.new.predict <- predict(umap.orig, out.lda.predict$topics)

umap.new.predict.long <- data.table(umap1 = umap.new.predict[, 1], umap2 = umap.new.predict[, 2], samp = rownames(umap.new.predict), batch = "Linneg") %>%
  mutate(umap2 = umap2 * -1)

ggplot(umap.new.predict.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

umap.merged <- bind_rows(umap.orig.long, umap.new.predict.long) %>%
  dplyr::rename(cell = samp)


ggplot(umap.merged, aes(x = umap1, y = umap2, color = batch)) + geom_point(alpha = 0.7, size = 1.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("H3K4me3")


# Get raw data  -----------------------------------------------------------

inf.raw <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-Linneg/PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.2019-06-17.merge_with_B6.RData"
load(inf.raw, v=T)


# Load proportions --------------------------------------------------------

probs.lst.filt <- readRDS("/Users/yeung/data/scchic/public_data/proportions_vectors_Lara-Astiaso.H3K4me3.2019-09-15.rds")


# Add decoys --------------------------------------------------------------

terms.keep <- names(probs.lst.filt[[1]])
# names(terms.keep) <- terms.keep

p.avg <- rowMeans(as.data.frame(probs.lst.filt))
p.flat <- rep(1 / length(p.avg), length(p.avg))

probs.lst.filt$Avg <- p.avg
probs.lst.filt$Unif <- p.flat


# Set up likelihoods ------------------------------------------------------


count.filt <- count.dat$counts
count.filt <- count.filt[terms.keep, ]

all.cells <- colnames(count.dat$counts)
names(all.cells) <- colnames(count.dat$counts)
cell.names <- names(all.cells)

cell.counts <- colSums(count.dat$counts) / 5
cell.counts.downsamp <- colSums(count.filt) / 5

LL.ctype.lst <- lapply(all.cells, function(cell.name){
  cell.vec <- count.filt[terms.keep, cell.name]
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
  })
})
# calculate probability of model given data 
p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)

# summaize
LL.dat <- lapply(cell.names, function(cname){
  LL.vec <- LL.ctype.lst[[cname]]
  p.vec <- p.ctype.lst[[cname]]
  cell.count = cell.counts[[cname]]
  cell.count.downsamp = cell.counts.downsamp[[cname]]
  if (all(is.infinite(LL.vec))){
    LL.max <- NA
    p.max <- NA
    best.ctype <- NA
  } else {
    LL.max <- max(LL.vec)
    p.max <- max(p.vec)
    best.ctype <- names(which.max(LL.vec))
  }
  dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, cell.count.downsamp = cell.count.downsamp)
  return(dat.tmp) 
}) %>%
  bind_rows()

LL.sum <- LL.dat %>%
  group_by(ctype.pred) %>%
  summarise(ncell = length(ctype.pred))

LL.dat.merge <- left_join(umap.merged, LL.dat)

# Plot again --------------------------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

LL.dat.merge.pretty <- LL.dat.merge %>%
  rowwise() %>%
  mutate(ctype.pred = ifelse(ctype.pred == "EryA" | ctype.pred == "EryB", "EryAorB", ctype.pred),
         ctype.pred = ifelse(ctype.pred == "CD4" | ctype.pred == "CD8", "CD4or8", ctype.pred))

# use cutoff to be confident on celltype prediction 


m.pretty <- ggplot(LL.dat.merge.pretty,
                   aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point(alpha=0.75, size = 3) + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark, "Probabilistic celltype inference using multinomial likelihood, \n with Ido Amit ChIP-seq data as underlying genomic region proportions") + 
  facet_wrap(~batch)
print(m.pretty)

# find an MF and Mono, whicih they different from Granu and B?
# take a low probaibility MF

# take all mono and arrange by rank?
jctype <- "Mono"
jctype <- "MF"
jctype <- "GN"
jctype <- "B"
jsub.not <- subset(LL.dat.merge, ctype.pred != jctype) %>% arrange(p.max) %>%
  mutate(rnk = 0,
         p.adj = NA)
jsub <- subset(LL.dat.merge, ctype.pred == jctype) %>% arrange(p.max) %>%
  mutate(rnk = seq(length(cell)),
         p.adj = ifelse(p.max > -0.1, p.max, NA))
jsub.merged <- bind_rows(jsub, jsub.not)
ggplot(jsub.merged, aes(x = umap1, y = umap2, color = rnk)) + 
  geom_point(alpha = 0.9) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("H3K4me3") + 
  scale_color_viridis_c()
ggplot(jsub.merged, aes(x = umap1, y = umap2, color = p.adj)) + 
  geom_point(alpha = 0.9) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("H3K4me3") + 
  scale_color_viridis_c(na.value = "gray85")
plot(jsub$rnk, jsub$p.max)

jsub <- subset(LL.dat.merge, ctype.pred == "Mono") %>% arrange(desc(p.max))
jcell <- jsub$cell[[1]]
Lvec <- sort(LL.ctype.lst[[jcell]], decreasing = TRUE)
Pvec <- sort(p.ctype.lst[[jcell]], decreasing = TRUE)
plot(Lvec)
text(Lvec, labels = names(Lvec))

plot(Pvec)
text(Pvec, labels = names(Lvec))

# plot on umap
ggplot(umap.merged %>% rowwise() %>% mutate(is.cell = cell == jcell), aes(x = umap1, y = umap2, color = is.cell, size = is.cell)) + 
  geom_point(alpha = 0.7) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("H3K4me3")


# plot probability?
m.pretty.prob <- ggplot(LL.dat.merge.pretty %>% filter(batch == "Ctrl"),
                   aes(x = umap1, y = umap2, color = p.max)) + geom_point(alpha=0.75, size = 3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark, "probabilistic celltype inference using multinomial likelihood, \n with ido amit chip-seq data as underlying genomic region proportions") + 
  facet_grid(batch~ctype.pred) +
  scale_color_viridis_c()
print(m.pretty.prob)

# plot probability?
m.pretty.prob <- ggplot(ll.dat.merge.pretty,
                   aes(x = umap1, y = umap2, color = p.max)) + geom_point(alpha=0.75, size = 3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(jmark, "probabilistic celltype inference using multinomial likelihood, \n with ido amit chip-seq data as underlying genomic region proportions") + 
  facet_grid(batch~ctype.pred) +
  scale_color_viridis_c()
print(m.pretty.prob)

# calculate enrichment of each celltype 
celltype.fc <- LL.dat.merge %>%
  group_by(batch, ctype.pred) %>%
  summarise(cell.count = length(cell)) %>%
  group_by(batch) %>%
  mutate(cell.frac = cell.count / sum(cell.count))

jctypes <- unique(celltype.fc$ctype.pred)
names(jctypes) <- jctypes
hyp.test.lst <- lapply(jctypes, DoFishersTestCtype, celltype.fc)

# plot odds ratio and p-value
hyp.test.dat <- lapply(jctypes, function(jctype){
  hyp.test <- hyp.test.lst[[jctype]]
  if (is.null(hyp.test)){
    warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
    return(NULL)
  }
  data.frame(OR = hyp.test$estimate, pval = hyp.test$p.value, ctype = jctype, stringsAsFactors = FALSE)
}) %>%
  bind_rows() %>%
  arrange(desc(OR)) %>%
  mutate(ctype = as.factor(ctype)) %>%
  dplyr::rename(ctype.pred = ctype)
  
m.enrichment <- ggplot(hyp.test.dat, aes(x = log10(OR), y = -log10(pval), label = ctype.pred)) + geom_point(size = 2.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_text_repel(size = 4) + 
  ggtitle("Fisher's exact test to quantify cell-type enrichment")
print(m.enrichment)

# sort by enrichment
celltype.fc.merge <- left_join(celltype.fc, hyp.test.dat)
m.barplot <- ggplot(celltype.fc.merge, aes(x = forcats::fct_reorder(.f = ctype.pred, .x = OR, .desc = TRUE), y = cell.frac, group = batch, fill = batch)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Celltype counts lineage-neg versus control, ordered by decreasing odds ratio")
print(m.barplot)

