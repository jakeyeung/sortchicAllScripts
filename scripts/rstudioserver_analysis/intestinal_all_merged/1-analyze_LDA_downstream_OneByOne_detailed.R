# Jake Yeung
# Date of Creation: 2020-01-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-analyze_LDA_downstream_OneByOne_detailed.R
# Look at things one by one because there is so much data

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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123



# Load public data  -------------------------------------------------------

inf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"

inf.pseudobulk <- "/home/jyeung/hpc/intestinal_scchic/public_data/pseudobulk_data/pseudobulk_Haber_et_al_intestinal_celltypes_2019-12-30.VST_noQuantNorm.rds"
dat.pseudobulk <- readRDS(inf.pseudobulk)

# do quant normalization?
library(preprocessCore)
dat.pseudobulk.norm <- normalize.quantiles(dat.pseudobulk, copy = TRUE)
rownames(dat.pseudobulk.norm) <- rownames(dat.pseudobulk)
colnames(dat.pseudobulk.norm) <- colnames(dat.pseudobulk)

boxplot(dat.pseudobulk)
boxplot(dat.pseudobulk.norm)

# dat.norm.long <- tidyr::gather(data.frame(gene = rownames(dat.pseudobulk.norm), dat.pseudobulk.norm, stringsAsFactors = FALSE), key = "celltype", value = "exprs", -gene) %>%
dat.norm.long <- tidyr::gather(data.frame(gene = rownames(dat.pseudobulk), dat.pseudobulk, stringsAsFactors = FALSE), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

# randomly select 150 genes and look at zscore distribution
set.seed(0)
jgenes.random <- sample(unique(dat.norm.long$gene), size = 150)
ggplot(subset(dat.norm.long, gene %in% jgenes.random), aes(x = celltype, y = zscore)) + geom_boxplot() + geom_jitter(width = 0.2)


# Load LDA  ---------------------------------------------------------------

jmark <- "k4me1"
jprefix <- "Scraped.AllMerged"
# jprefix <- "Scraped.Unenriched"
jbin <- FALSE
jdate <- "2019-12-22"

inf <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.", jprefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.", jdate, ".K-30.binarize.", jbin, "/ldaOut.mat.", jprefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.", jdate, ".K-30.Robj")
if (!file.exists(inf)){
  jdate <- "2019-12-23"
  print(paste("Using date:", jdate))
  inf <- paste0("/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22/lda_outputs.mat.", jprefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.", jdate, ".K-30.binarize.", jbin, "/ldaOut.mat.", jprefix, ".", jmark, ".TAcutoff_0.5.countscutoff_500_1000.binfilt_cellfilt.", jdate, ".K-30.Robj")
}
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

tm.result <- posterior(out.lda)
colnames(tm.result$topics) <- paste("Topic", colnames(tm.result$topics), sep = "_")

# annotate bins
annot.out <- AnnotateBins(terms.mat = tm.result$terms, inf.tss = inf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")

terms.sum <- annot.out$terms.filt %>%
  group_by(gene) %>%
  dplyr::filter(rnk == min(rnk))

term2gene <- hash(terms.sum$term, terms.sum$gene)
gene2term <- hash(terms.sum$gene, terms.sum$term)

terms.filt <- annot.out$terms.filt %>%
  mutate(topic = paste("Topic_", topic, sep = ""))

topics.mat <- tm.result$topics

umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell),
         prefix = paste(strsplit(cell, "-")[[1]][1:3], collapse = "-"))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = prefix)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~prefix)

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi)

# plot variance
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)
dat.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + scale_color_viridis_c(direction = -1)


# what are the celltypes in the islands????

# get different topics
topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)

dat.topics <- left_join(dat.umap.long, data.frame(cell = rownames(tm.result$topics), tm.result$topics, stringsAsFactors = FALSE))

jtop <- topics.sum$topic[[1]]
jtop <- topics.sum$topic[[2]]
jtop <- topics.sum$topic[[3]]




# Do all topics -----------------------------------------------------------

keeptop <- 100

jtop <- topics.sum$topic[[4]]


m.umap <- PlotXYWithColor(dat.topics, xvar = "umap1", yvar = "umap2", cname = jtop)


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



top.genes <- subset(terms.filt, topic == jtop & rnk <= keeptop)$gene
assertthat::assert_that(length(top.genes) > 0)

jsub <- subset(dat.norm.long, gene %in% top.genes)
jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% 
  summarise(zscore = median(zscore)) %>% 
  arrange(desc(zscore)) %>% dplyr::select(celltype)
jlevels <- as.character(jsub.sorted.summarised$celltype)
jsub$celltype <- factor(jsub$celltype, levels = jlevels)

m.exprs <- ggplot(jsub,
                  aes(x = celltype , y = zscore)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin() +
  geom_jitter(width = 0.1, size = 0.5) +
  # geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))

print(m.umap)
print(m.top)
print(m.exprs)



# Plot marker genes  ------------------------------------------------------

jdist <- "50000"
inf.mat.tss <- paste0("/home/jyeung/hpc/intestinal_scchic/raw_data/HVG-intestines.tagged_bams.all_merged/countTables_TSS/HVG_Scraped_Bl6_k4me1_2019-12-20.TSS_", jdist, ".countTable.csv")
assertthat::assert_that(file.exists(inf.mat.tss))
count.mat.tss <- ReadMatTSSFormat(inf.mat.tss)

count.mat.tss.long <- CollapseRowsByGene(count.mat.tss, as.long = TRUE)

colnames(count.mat.tss.long) <- c("gene", "cell", "count")


# normalize by total
count.mat.tss.long <- count.mat.tss.long %>%
  group_by(cell) %>%
  mutate("frac.umi" = (count / sum(count)),
         "log2rpm" = log2(frac.umi * 10^6 + 1))

jgenes <- c("Alpi", "Lyz1", "Chgb", "Clca3", "Lgr5", "Msi2", "Hoxb5", "Cerkl", "Ankrd11", "Gimap6", "Prkch", "Mbnl1", 
            "Kalrn", "Apob", "Jarid2", "Arhgap42", "Gsk3b",  # entrocytes?
            "Defa23", "Defa17", "Defa3", "Defa22",  # Paneth
            "Muc2", "Clca3", "S100a6", "Scin",  # Goblet??
            "Foxa2")  

# choose a few genes and integrate into matrix
dat.exprs.sub <- spread(subset(count.mat.tss.long, gene %in% jgenes, select = c(cell, gene, log2rpm)), key = "gene", value = "log2rpm")

dat.exprs.merge <- left_join(dat.umap.long, dat.exprs.sub)

jgene <- "Cerkl"
PlotXYWithColor(dat.exprs.merge, xvar = "umap1", yvar = "umap2", cname = jgene) + scale_color_viridis_c()


# Show louvains  ----------------------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)

# try a Poisson regression model 

# try enteroendocrine?
jlouv <- 1

jlouvs <- unique(dat.umap.long$louvain)

jstart <- Sys.time()
for (jlouv in jlouvs){
  print(paste("Fits for:", jlouv))
  
  outmain <- "/home/jyeung/data/from_rstudioserver/scchic/intestinal_analysis"
  outpdf <- file.path(outmain, paste0("louv_plot_", jlouv, ".pdf"))
  
  jsub.louv <- dat.umap.long %>% rowwise() %>% mutate(louvain = ifelse(louvain == jlouv, "InGroup", "OutGroup"))
  m.louv <- ggplot(jsub, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + ggtitle(paste("Louvain:", jlouv))
  
  pdf(outpdf, useDingbats = FALSE)
    print(m.louv)
  dev.off()
  
  outf <- file.path(outmain, paste0("glm_fits_louv_", jlouv, ".rds"))
  if (file.exists(outf)){
    next
  }
  
  # do glm genome wide
  count.mat.tss.long.filt.annot <- left_join(count.mat.tss.long, subset(jsub.louv, select = c(cell, louvain))) %>%
    filter(!is.na(louvain))
  # fit a gene
  # do genome wide
  system.time(
    jfits <- count.mat.tss.long.filt.annot %>%
      group_by(gene) %>%
      do(fitGene = glm(formula = count ~ louvain, family = "poisson", data = .))
  )
  saveRDS(jfits, file = outf)
}
print(Sys.time() - jstart)


# check cell not in group
cells.not.in <- colnames(count.mat.tss)[!(colnames(count.mat.tss) %in% dat.umap.long$cell)]
subset(count.mat.tss.long.annot, cell == cells.not.in[[1]])


