# Jake Yeung
# Date of Creation: 2019-10-29
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/3-celltype_analysis_enrichment_with_transcriptome_redo.R
# Redo to repdocue


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(hash)
library(igraph)
library(umap)
library(Seurat)
library(scchicFuncs)

library(preprocessCore)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DESeq2)


# Functions ---------------------------------------------------------------



# Functions ---------------------------------------------------------------

DoFishersTestCtype <- function(jctype, celltype.fc){
  conting.tab <- spread(celltype.fc %>%
                          mutate(ctype.pred = ifelse(ctype.pred == jctype, jctype, "zNot")) %>%
                          dplyr::select(is.stem, ctype.pred, cell.count) %>%
                          group_by(ctype.pred, is.stem) %>%
                          summarise(cell.count = sum(cell.count)),
                        key = is.stem, value = cell.count) %>%
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


DoFishersTestCtypeStringent <- function(jctype, celltype.fc){
  conting.tab <- spread(celltype.fc %>%
                          mutate(ctype.stringent = ifelse(ctype.stringent == jctype, jctype, "zNot")) %>%
                          dplyr::select(is.stem, ctype.stringent, cell.count) %>%
                          group_by(ctype.stringent, is.stem) %>%
                          summarise(cell.count = sum(cell.count)),
                        key = is.stem, value = cell.count) %>%
    as.data.frame() %>%
    ungroup()
  if (any(is.na(conting.tab))){
    warning(paste0("celltype:", jctype, ". NAs found, probably no counts in one of the tables, returning NULL"))
    return(NULL)
  }
  # rownames
  rownames(conting.tab) <- conting.tab$ctype.stringent
  conting.tab$ctype.stringent <- NULL
  
  # make odds ratio interpretable:( zNot_treat / ctype_treat ) / ( zNot_ctrl / ctype_ctrl )
  conting.tab <- t(conting.tab)
  # swap rows
  conting.tab <- conting.tab[c(2, 1), ]
  print(conting.tab)
  hyp.test <- fisher.test(x = conting.tab)
  return(hyp.test)
}


# Set up tables -----------------------------------------------------------


# plot out
jsettings <- umap.defaults
jsettings$n_neighbors <- 25
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.Robj"
load(inf.lda, v=T)

out.lda <- out.lda[[3]]

umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
dat.umap.long.lda <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

dat.umap.long.lda$is.stem <- sapply(dat.umap.long.lda$cell, function(x) grepl("stem-cell", x))

ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2)) + geom_point() + 
  facet_wrap(~is.stem) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.long.lda <- DoLouvain(topics.mat = posterior(out.lda)$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.lda)

ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  facet_wrap(~is.stem) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


inf.proj <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/projections/H3K4me1_matsMergedNonenriched.matsMergedEnriched.RData"
load(inf.proj, v=T)

dat.pred <- predict(umap.out, out.lda.predict$topics)

dat.pred.long <- data.frame(cell = rownames(dat.pred), umap1 = dat.pred[, 1], umap2 = dat.pred[, 2], stringsAsFactors = FALSE) %>%
  mutate(is.stem = TRUE)

dat.umap.pred.merged <- bind_rows(subset(dat.umap.long.lda, select = -louvain), dat.pred.long)

ggplot(dat.umap.pred.merged, aes(x = umap1, y = umap2)) + geom_point() + 
  facet_wrap(~is.stem) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Get raw counts ----------------------------------------------------------

inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedAll_2019-09-29.CountThres0.K-30_35_50.OutObjs.RData"
load(inf.lda, v=T)


# Layer on K4me1 data -----------------------------------------------------

# load public data 


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  group_by(Gene_Name, CellType) %>%
  summarise(FPKM = sum(FPKM)) %>%
  rowwise() %>%
  mutate(logFPKM = log2(FPKM + 1))

# normalize across samples?
ggplot(dat.long, aes(x = CellType, y = logFPKM)) + geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dat.mat <- tidyr::spread(dat.long %>%
                           ungroup() %>%
                           # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>%
                           mutate(gene = Gene_Name) %>%
                           dplyr::select(gene, CellType, logFPKM),
                         key = CellType, value = logFPKM)  %>%
  as.data.frame()
rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL

cnames.tmp <- colnames(dat.mat)
rnames.tmp <- rownames(dat.mat)
dat.mat <- preprocessCore::normalize.quantiles(as.matrix(dat.mat), copy = TRUE)  # strong normalization,
colnames(dat.mat) <- cnames.tmp
rownames(dat.mat) <- rnames.tmp

boxplot(dat.mat)

dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>%
  as.data.frame()
rownames(dat.norm.zscore.mat) <- dat.norm.zscore.mat$gene
dat.norm.zscore.mat$gene <- NULL

annot.out <- AnnotateBins(terms.mat = out.objs$tm.result$terms)

terms.sum <- annot.out$terms.filt %>%
  group_by(gene) %>%
  dplyr::filter(rnk == min(rnk))

term2gene <- hash(terms.sum$term, terms.sum$gene)
gene2term <- hash(terms.sum$gene, terms.sum$term)

genes.keep <- rownames(dat.mat)

top.nterms <- 150
terms.all <- unique(subset(terms.sum %>% arrange(desc(weight)), rnk <= top.nterms)$term)
genes.all <- unlist(sapply(terms.all, function(x) term2gene[[x]]))



# Do likelihoods ----------------------------------------------------------

# create likelihoods
jcutoff <- 1.8
plot(density(rowMeans(dat.mat)))
abline(v = jcutoff)

dat.mat.filt <- dat.mat[apply(dat.mat, 1, max) > jcutoff, ]

genes.keep <- intersect(genes.all, rownames(dat.mat.filt))
terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])

dat.mat.filt <- dat.mat.filt[genes.keep, ]



# redefine genes that we keep
genes.filt <- intersect(rownames(dat.mat.filt), genes.keep)

# dat.mat.filt <- dat.mat.filt[genes.filt, ]
# try using zscore
dat.mat.filt <- 2^dat.norm.zscore.mat[genes.filt, ]  # exp is worse than 2

# dat.mat.filt <- exp(dat.norm.zscore.mat[genes.filt, ])
# renormalize so no zeros??
# dat.mat.filt <- sweep(dat.mat.filt, MARGIN = 2, STATS = apply(dat.mat.filt, 2, min), FUN = "-")
# dat.mat.filt <- as.matrix(dat.mat.filt)


# handle zeros
zero.fill <- min(as.matrix(dat.mat.filt)[which(as.matrix(dat.mat.filt) != 0)])
dat.mat.filt[which(dat.mat.filt == 0)] <- zero.fill

# make likelihoods
probs.lst.filt <- as.list(as.data.frame(dat.mat.filt))
# name the list just to be safe
probs.lst.filt <- lapply(probs.lst.filt, function(x){
  names(x) <- rownames(dat.mat.filt)
  return(x)
}) 


# Prepare count mat -------------------------------------------------------

count.filt <- out.objs$count.mat[terms.keep, ]
rownames(count.filt) <- genes.keep

# Do fits -----------------------------------------------------------------

all.cells <- colnames(out.objs$count.mat)
names(all.cells) <- all.cells

LL.ctype.lst <- lapply(all.cells, function(cell.name){
  cell.vec <- count.filt[genes.filt, cell.name]
  # cell.vec <- cell.vec[which(cell.vec > 0)]
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    # jprob <- jprob[names(cell.vec)]
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    # return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
    # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
    # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
    # return(dmultinom(x = cell.vec, prob = sqrt(jprob), log = TRUE))
    return(dmultinom(x = cell.vec, prob = jprob^(1/4), log = TRUE))
  })
})

# Summarize fits ----------------------------------------------------------

# calculate probability of model given data 
p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)

cell.counts <- Matrix::colSums(out.objs$count.mat) / 5
cell.counts.downsamp <- Matrix::colSums(out.objs$count.mat) / 5

cell.names <- names(all.cells)

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
  dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, cell.count.downsamp = cell.count.downsamp, stringsAsFactors = FALSE)
  return(dat.tmp) 
}) %>%
  bind_rows()

# be stringent with the predictions 


LL.sum <- LL.dat %>%
  group_by(ctype.pred) %>%
  summarise(ncell = length(ctype.pred))

print(LL.sum)

# p.filt <- log(0.9)
p.filt <- log(0)
LL.dat <- LL.dat %>%
  rowwise() %>%
  mutate(ctype.stringent = ifelse(p.max >= p.filt, ctype.pred, NA))

LL.dat.merge <- left_join(dat.umap.pred.merged , LL.dat)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.umap.celltype <- ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = ctype.stringent)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey95") +
  theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.stringent) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



ggplot(LL.dat.merge, aes(x = umap1, y = umap2, color = LL.max/cell.size)) + geom_point() + 
  scale_color_viridis_c() + 
  facet_wrap(~ctype.stringent) + theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show enrihment before and after
m.enrichment.umap <- ggplot(subset(LL.dat.merge, !is.na(ctype.stringent)), aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~ctype.stringent)

print(m.enrichment.umap + ggtitle(paste("Pfilt", p.filt)))



# Get enrichment ----------------------------------------------------------

# calculate enrichment of each celltype
celltype.fc <- LL.dat.merge %>%
  filter(!is.na(ctype.stringent)) %>%
  group_by(is.stem, ctype.stringent) %>%
  summarise(cell.count = length(cell)) %>%
  group_by(is.stem) %>%
  mutate(cell.frac = cell.count / sum(cell.count))

jctypes <- unique(celltype.fc$ctype.stringent)
names(jctypes) <- jctypes
hyp.test.lst <- lapply(jctypes, DoFishersTestCtypeStringent, celltype.fc)

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
  dplyr::rename(ctype.stringent = ctype)

m.enrichment <- ggplot(hyp.test.dat, aes(x = log10(OR), y = -log10(pval), label = ctype.stringent)) + geom_point(size = 2.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text_repel(size = 4) +
  ggtitle("Fisher's exact test to quantify cell-type enrichment")
print(m.enrichment)

# sort by enrichment
celltype.fc.merge <- left_join(celltype.fc, hyp.test.dat)
m.barplot <- ggplot(celltype.fc.merge, aes(x = forcats::fct_reorder(.f = ctype.stringent, .x = OR, .desc = TRUE), y = cell.frac, group = is.stem, fill = is.stem)) + geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Celltype counts lineage-neg versus control, ordered by decreasing odds ratio", p.filt)
print(m.barplot)


# write outputs
outpdf <- paste0("/Users/yeung/data/scchic/pdfs/stemcell_analysis/celltype_enrichment_with_tx.", Sys.Date(), ".pdf")

pdf(outpdf, useDingbats = FALSE)

  print(m.umap.celltype)
  print(m.umap.celltype + facet_wrap(~is.stem))
  print(m.enrichment.umap)
  print(m.enrichment)
  print(m.barplot)

dev.off()

