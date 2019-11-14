# Jake Yeung
# Date of Creation: 2019-11-13
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_TSS_multinom_celltyping.R
# Multinomial celltyping

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

library(ggrepel)

source("/Users/yeung/projects/scchicFuncs/R/MultinomFunctions.R")


# Load data ---------------------------------------------------------------

jprefix <- "ZFWKM"

# jmarks <- c("H3K4me1", "H3K4me3")
# winsizes <- c(50000L, 100000L)

jmark <- "H3K4me3"
# winsize <- 50000L
winsize <- 100000L

inf.annot <- paste0("/Users/yeung/data/scchic/tables/gene_tss.winsize_", winsize, ".species_drerio.nochr.bed")

# init
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-", jprefix, "-", jmark, 
              ".winsize_", winsize, ".merged.K-30.binarize.FALSE/ldaOut.PZ-ChIC-", jprefix, "-", jmark, 
              ".winsize_", winsize, ".merged.K-30.Robj")
assertthat::assert_that(file.exists(inf))

x <- load(inf, v=T)

if (length(out.lda) > 1){
  out.lda <- out.lda[[1]]
} 

# prepare raw data 
genes.chic <- sapply(rownames(count.mat), function(x) strsplit(x, ";")[[1]][[2]])
# get count.filt, sum across same gene
count.mat.tmp <- count.mat
rownames(count.mat.tmp) <- genes.chic
count.mat.long <- melt(as.matrix(count.mat.tmp), value.name = "count") %>%
  group_by(Var1, Var2) %>%
  summarise(count = max(count))
  # summarise(count = sum(count))
rm(count.mat.tmp)


# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
inf.WKM <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"
dat.bulk <- readRDS(inf.WKM)
dat.bulk.mat <- dcast(subset(dat.bulk, select = c(gene, celltype, exprs)), gene ~ celltype, value.var = "exprs")
rownames(dat.bulk.mat) <- dat.bulk.mat$gene; dat.bulk.mat$gene <- NULL

zscore.cutoff <- 1
zscore.cutoff <- "top500"
# set up reference data

if (is.character(zscore.cutoff)){
dat.bulk.keep <- dat.bulk %>%
  group_by(celltype) %>%
  mutate(rnk = rank(-zscore)) %>%
  filter(rnk < 500)
} else {
  dat.bulk.keep <- dat.bulk %>%
    group_by(gene) %>%
    filter(max(abs(zscore)) > zscore.cutoff)
}

ref.genes.keep <- intersect(as.character(dat.bulk.keep$gene), genes.chic)

print(length(ref.genes.keep))


count.filt <- dcast(subset(count.mat.long, Var1 %in% ref.genes.keep), Var1 ~ Var2, value.var = "count")
rownames(count.filt) <- count.filt$Var1; count.filt$Var1 <- NULL
count.filt <- as.matrix(count.filt)[ref.genes.keep, ]  # keep same order!

# Analyze LDA -------------------------------------------------------------


# add gene name to the coordinates (got lost in mat to sparse mat pipeline)
topics.mat <- posterior(out.lda)$topics
terms.mat <- posterior(out.lda)$terms

colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))
rownames(terms.mat) <- paste0("topic_", rownames(terms.mat))

print(head(out.lda@terms))

# # if no gene names in rownames, then add them 
# annots.hash <- GetGeneAnnotsHash(inf.annot)
# count.mat <- AddGeneNameToRows(count.mat, annots.hash)

# out.lda@terms <- sapply(out.lda@terms, function(x) annots.hash)


# colnames(terms.mat) <- sapply(colnames(terms.mat), function(x) annots.hash[[x]])

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jsettings.1d <- jsettings; jsettings.1d$n_components <- 1

umap.out <- umap(topics.mat, config = jsettings)
umap.out.1d <- umap(topics.mat, config = jsettings.1d)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
dat.umap.long.1d <- data.frame(cell = rownames(umap.out$layout), umap1.1d = umap.out$layout[, 1])

m.umap.blank <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("LDA with windows around TSS. Winsize:", winsize))

# get variance??
dat.impute.log <- log2(t(topics.mat %*% terms.mat))
rownames(dat.impute.log) <- gsub(";", "_", rownames(dat.impute.log))

# intrachromosomal variance doesnt make sense when doing TSS, do genome wide
# jchromos <- paste("chr", seq(25), sep = "")
jchromos <- c("")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.long.varmerge <- left_join(dat.umap.long, dat.var)
dat.umap.long.varmerge <- left_join(dat.umap.long.varmerge, dat.umap.long.1d)

PlotXYWithColor(dat.umap.long.varmerge, xvar = "umap1", yvar = "umap2", cname = "cell.var.within.sum.norm")

ggplot(dat.umap.long.varmerge, aes(x = cell.var.within.sum.norm, y = umap1.1d)) + geom_point() 


print(m.umap.blank)

m.umap.var <- ggplot(dat.umap.long.varmerge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()
print(m.umap.var)


# Find celltypes ----------------------------------------------------------

# get topics
topics.sum <- OrderTopicsByEntropy(posterior(out.lda), jquantile = 0.99) %>%
  rowwise() %>%
  mutate(topic = gsub("^X", "topic_", topic))

# show a topic

dat.umap.long.merge <- left_join(dat.umap.long, data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE))
colnames(dat.umap.long.merge) <- gsub("^X", "topic_", colnames(dat.umap.long.merge))

jtop <- topics.sum$topic[[1]]

m.umap.top <- PlotXYWithColor(dat.umap.long.merge, xvar = "umap1", yvar = "umap2", cname = jtop)

# what are the genes? 

# Plot Genes Across Topic

# make terms mat once
terms.mat.long <- data.table::melt(terms.mat)
colnames(terms.mat.long) <- c("topic", "term", "weight")
terms.mat.long$term <- as.character(terms.mat.long$term)
terms.mat.long <- terms.mat.long %>%
  rowwise() %>%
  mutate(gene = strsplit(term, ";")[[1]][[2]]) %>%
  group_by(topic) %>%
  mutate(rnk = rank(-weight)) %>%
  arrange(desc(weight))


# set up probs.lst.filt

probs.lst.filt <- SetUpProbs(dat.bulk.mat[ref.genes.keep, ])

all.cells <- colnames(count.filt)
names(all.cells) <- all.cells

LL.ctype.lst <- FitMultinoms(count.filt, all.cells, probs.lst.filt, exppower = 0.25)



# Summarize fits ----------------------------------------------------------

LL.dat <- SummarizeMultinomFits(LL.ctype.lst, count.filt)

# add to umap 

dat.umap.long.merge <- left_join(dat.umap.long, LL.dat, by = "cell")

jtitle <- paste(paste0("Mark_", jmark), paste0("winsize_", winsize), paste0("Zcutoff_", zscore.cutoff), sep = ".")
m.celltypes <- ggplot(dat.umap.long.merge %>% filter(p.max > log(0.9)), aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~ctype.pred) + ggtitle(jtitle) + 
  scale_color_viridis_c()

outdir <- "/Users/yeung/data/scchic/pdfs/zebrafish/celltyping"
pdf(file.path(outdir, paste0(jtitle, ".pdf")), useDingbats = FALSE)
  print(m.celltypes)
dev.off()

print(m.celltypes)
