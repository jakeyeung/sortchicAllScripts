# Jake Yeung
# Date of Creation: 2019-12-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/3-celltyping_on_LDA.R
# 



rm(list=ls())

jstart <- Sys.time()

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

library(ggrepel)

library(here)

library(topicmodels)

setwd(here())

# Load data  --------------------------------------------------------------

# jmark <- "H3K27me3"

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jsuff <- "BM"
jbin <- "FALSE"
for (jmark in jmarks){
  

# outmain <- "/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs"
outmain <- "/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_NoSlideWin"
dir.create(outmain)

infmain <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE")
assertthat::assert_that(dir.exists(infmain))

inf.bulkdat <- "/home/jyeung/hpc/scChiC/public_data/E-MTAB-3079-query-results.fpkms.tsv"
assertthat::assert_that(file.exists(inf.bulkdat))

inf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"

keeptop <- 150


outpdfname <- paste0(jmark, ".", jsuff, ".bin_", jbin, ".celltype_using_topics.pdf")
outpdf <- file.path(outmain, outpdfname)
if (file.exists(outpdf)){
  print(paste("Outpdf exists, skipping for", jmark))
  next
}

inf <- file.path(infmain, paste0("ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj"))

load(inf, v=T)

if (is.na(all(count.mat.orig))){
  count.mat.orig <- count.mat
}

tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)  %>%
  rowwise() 

dat.umap.long <- dat.umap.long %>%
  mutate(experi = ClipLast(cell, jsep = "_"))

print(unique(dat.umap.long$experi))

assertthat::assert_that(length(unique(dat.umap.long$experi)) < 20)  # need to double check experi

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

dat.impute.log <- log2(t(topics.mat %*% terms.mat))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

# Load bulk ---------------------------------------------------------------

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

annot.out <- AnnotateBins2(terms.mat = terms.mat, inf.tss = inf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")

annot.out$terms.annot <- annot.out$terms.annot %>%
  filter(!is.na(termgene)) %>%
  mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
  group_by(gene)  # dont filter use

terms.sum <- annot.out$terms.annot %>%
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

dat.mat.filt <- dat.mat.filt[genes.filt, ]

probs.lst.filt <- SetUpProbs(dat.mat.filt, norm.vec = TRUE)

# Prepare count mat -------------------------------------------------------

count.filt <- count.mat.orig[terms.keep, ]
rownames(count.filt) <- genes.keep


# Run fits ----------------------------------------------------------------


all.cells <- colnames(count.filt)
names(all.cells) <- all.cells

LL.ctype.lst <- FitMultinoms(count.filt, all.cells, probs.lst.filt, exppower = 0.25)


# Summarize ---------------------------------------------------------------

LL.dat <- SummarizeMultinomFits(LL.ctype.lst, count.filt, all.cells) %>%
  mutate(is.stem = grepl("Linneg", cell))

m.celltypes.bar <- ggplot(LL.dat %>% group_by(is.stem, ctype.pred) %>% summarise(count = length(ctype.pred)) %>% group_by(is.stem) %>% mutate(total = sum(count)),
                          aes(x = ctype.pred, y = count / total, group = is.stem, fill = is.stem)) +
  geom_bar(position = "dodge", stat = "identity")  +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(m.celltypes.bar)

# Do UMAP  ----------------------------------------------------------------

dat.merged <- left_join(dat.umap.long, LL.dat)
dat.merged <- left_join(dat.merged, dat.var)

dat.merged <- dat.merged %>%
  mutate(ctype.pred.stringent = ifelse(p.max > log(0.99), ctype.pred, NA))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.merged, aes(x = umap1, y = umap2, color = experi)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette)

ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette)

ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred.stringent)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value="gray80") + facet_wrap(~ctype.pred)

# show predicted celltypes have low variance
ggplot(dat.merged, aes(x = cell.var.within.sum.norm)) + geom_density(alpha = 0.5, fill = 'blue')  + 
  theme_bw() + theme(aspect.ratio=0.1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~ctype.pred, ncol = 1)

# calculate average
dat.var.avg <- subset(dat.merged, !is.na(ctype.pred.stringent)) %>%
  group_by(ctype.pred.stringent) %>%
  summarise(var.mean = mean(cell.var.within.sum.norm))

ggplot(dat.var.avg, aes(x = forcats::fct_reorder(.f = ctype.pred.stringent, .x = var.mean, .fun = median, .desc = TRUE), y = var.mean)) + geom_bar(stat = "identity") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.merged, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)

ggplot(dat.merged, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~ctype.pred.stringent)



# Do celltyping based on topics  ------------------------------------------

tm.result.tmp <- tm.result
colnames(tm.result.tmp$topics) <- paste("Topic_", colnames(tm.result.tmp$topics), sep = "")

topics.sum <- OrderTopicsByEntropy(tm.result.tmp, jquantile = 0.99)

# add to dat
dat.merged.topics <- left_join(dat.umap.long, data.frame(cell = rownames(tm.result.tmp$topics), tm.result.tmp$topics))

terms.filt <- annot.out$terms.annot %>%
  mutate(topic = paste("Topic_", topic, sep = ""))

jtop <- "Topic_25"


pdf(outpdf, useDingbats = FALSE)
m.umap.experi <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.umap.experi2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

m.umap.var <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~experi)
m.umap.pred <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette)

m.umap.pred2 <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred.stringent)) + geom_point(alpha = 1, size = 5)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") +
  scale_color_manual(values = cbPalette)

m.umap.stringent <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred.stringent)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value="gray80") + facet_wrap(~ctype.pred)

# show predicted celltypes have low variance
m.umap.ctype.var <- ggplot(dat.merged, aes(x = cell.var.within.sum.norm)) + geom_density(alpha = 0.5, fill = 'blue')  + 
  theme_bw() + theme(aspect.ratio=0.1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~ctype.pred, ncol = 1)


print(m.umap.experi)
print(m.umap.var)
print(m.umap.pred)
print(m.umap.stringent)
print(m.umap.ctype.var)

for (jtop in topics.sum$topic){
  
  print(jtop)
  m.umap <- PlotXYWithColor(dat.merged.topics, xvar = "umap1", yvar = "umap2", cname = jtop)
  
  top.genes <- subset(terms.filt, topic == jtop & rnk <= keeptop)$gene
  assertthat::assert_that(length(top.genes) > 0)
  
  jsub <- subset(dat.norm.long, gene %in% top.genes)
  jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(celltype)
  jlevels <- as.character(jsub.sorted.summarised$celltype)
  jsub$celltype <- factor(jsub$celltype, levels = jlevels)
  m.exprs <- ggplot(jsub,
                    aes(x = celltype , y = zscore)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_violin() +
    geom_jitter(width = 0.1, size = 0.5) +
    # geom_line() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
    ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
  # plot top 150 genes?
  jsub.terms <- subset(terms.filt, topic == jtop & rnk < keeptop) %>%
    ungroup() %>%
    mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
  m.top <- jsub.terms %>%
    # mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
    ggplot(aes(x = term, y = log10(weight), label = gene)) +
    geom_point(size = 0.25) +
    theme_bw(8) +
    geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(paste("Top peak weights for:", jtop))
  
  # plot everything
  
  print(m.umap)
  print(m.exprs)
  print(m.top)
}
dev.off()

}


print(Sys.time() - jstart)
