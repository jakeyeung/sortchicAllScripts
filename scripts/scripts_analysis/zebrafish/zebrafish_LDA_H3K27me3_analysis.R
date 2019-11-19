# Jake Yeung
# Date of Creation: 2019-11-14
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_H3K27me3_analysis.R
# from # File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_merged_analysis.R
# merged analysis

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)
library(DESeq2)
library(preprocessCore)
library(ggrepel)
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)


# Load dat ----------------------------------------------------------------

jmark <- "H3K27me3"
jbin <- "FALSE"

# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.PZ-ChIC-ZFWKM-", jmark, ".2019-11-13.K-30_35_50.binarize.", jbin, "/ldaOut.PZ-ChIC-ZFWKM-", jmark, ".2019-11-13.K-30_35_50.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# Plot data ---------------------------------------------------------------

kchoose <- 30
kchoose.i <- which(sapply(out.lda, function(x) x@k) == kchoose)
topics.mat <- posterior(out.lda[[kchoose.i]])$topics
colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("CD41plus", cell))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)


# Do variance calculation -------------------------------------------------

jfac <- 10^6
jpseudo <- 0
dat.impute.log <- log2(t(posterior(out.lda[[1]])$topics %*% posterior(out.lda[[1]])$terms) * jfac + jpseudo)
jchromos.num <- seq(25)
jchromos <- paste("chr", jchromos.num, sep = "")

cells.var.chromo.merged <- CalculateVarAll(dat.impute.log, jchromos)


# Plot with variance  -----------------------------------------------------

dat.umap.long.merge <- left_join(dat.umap.long, cells.var.chromo.merged)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  scale_colour_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark) + 
  facet_wrap(~is.stem)

# do histogram?
ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm)) + geom_histogram() + facet_wrap(~is.stem, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark) + xlab("Intrachromosomal Variance")



# Get regions ------------------------------------------------------------

tm.result <- posterior(out.lda[[kchoose.i]])

regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd),
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
# remove weird chromo names
jchromos <- paste("chr", seq(25), sep = "")
regions <- subset(regions, seqnames %in% jchromos)

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))

regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                TxDb=TxDb.Drerio.UCSC.danRer11.refGene,
                                                annoDb='org.Dr.eg.db'))
regions.annotated$region_coord <- names(regions.range)
top.thres <- 0.99
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
}) %>%
  unlist() %>%
  unique()

regions.annotated.filt <- subset(regions.annotated, distanceToTSS < 1000)

bins.filt <- rownames(regions)

terms.dat <- data.frame(term = rownames(t(tm.result$terms)), t(tm.result$terms), stringsAsFactors = FALSE) %>%
  gather(., key = topic, value = weight, -term) %>%
  mutate(topic = gsub("^X", "topic_", topic)) %>%
  dplyr::filter(term %in% bins.filt) %>%
  group_by(topic) %>%
  arrange(desc(weight)) %>%
  mutate(jrank = seq(length(weight))) %>%
  mutate(term = as.factor(term)) %>%
  left_join(., subset(regions.annotated.filt, select = c(region_coord, SYMBOL)) %>% dplyr::rename(term = region_coord, gene = SYMBOL))



# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
inf.WKM <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"
dat.bulk <- readRDS(inf.WKM)
dat.bulk.mat <- dcast(subset(dat.bulk, select = c(gene, celltype, exprs)), gene ~ celltype, value.var = "exprs")
rownames(dat.bulk.mat) <- dat.bulk.mat$gene; dat.bulk.mat$gene <- NULL

dat.mat.filt.long <- tidyr::gather(data.frame(gene = rownames(dat.bulk.mat), dat.bulk.mat), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  ungroup()

# Get topics  -------------------------------------------------------------



# Look at topics and associate with celltypes -----------------------------

# order topics

topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
topics.sum$topic <- gsub("^X", "topic_", topics.sum$topic)

# show an island
topics.mat.tmp <- as.data.frame(topics.mat)
topics.mat.tmp$cell <- rownames(topics.mat)

topics.mat.tmp <- left_join(dat.umap.long, topics.mat.tmp)

keepns <- c(150)
# keepn <- 25

for (keepn in keepns){
  print(paste("keepn", keepn))
  outpdf <- paste0("/Users/yeung/data/scchic/pdfs/zebrafish/ZF_", jmark, "_topics_celltyping.keepn_", keepn, ".binarize_", jbin, ".kchoose_", kchoose, ".pdf")
  pdf(outpdf, useDingbats = FALSE)
  for (i in seq(nrow(topics.sum))){
    jtop <- topics.sum$topic[[i]]
    print(jtop)
    # plot UMAP
    m.umap <- PlotXYWithColor(topics.mat.tmp, xvar = "umap1", yvar = "umap2", cname = jtop, cont.color = TRUE, jtitle = jtop) + scale_color_viridis_c()
    
    jsub.terms <- subset(terms.dat, topic == jtop) %>%
      dplyr::filter(jrank <= keepn)
    
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
      ggtitle(paste("Top peak weights for:", jtop))
    
    # check gene expression across genes
    gfilt <- unique(jsub.terms$gene)
    m.ctype <- subset(dat.mat.filt.long, gene %in% gfilt) %>%
      mutate(celltype = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median)) %>%
      ggplot(., aes(x = celltype, y = zscore)) + geom_boxplot() + 
      geom_point() + 
      theme_bw(24) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ggtitle(jtop)
    
    print(m.umap)
    print(m.top)
    print(m.ctype)
  }
  dev.off()
}


