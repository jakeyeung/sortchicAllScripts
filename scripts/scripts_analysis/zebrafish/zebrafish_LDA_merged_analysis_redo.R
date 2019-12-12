# Jake Yeung
# Date of Creation: 2019-11-24
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_merged_analysis_redo.R
# 

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

library(xlsx)

library(DESeq2)

library(preprocessCore)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)

library(ggrepel)

# Load dat ----------------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/pdfs/zebrafish.redo"

jmark <- "H3K9me3"
jbin <- "TRUE"

# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_Merged.", jmark, ".2019-11-22.K-50.binarize.", jbin, "/ldaOut.PZ-ZF-All_Merged.", jmark, ".2019-11-22.K-50.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

# Plot data ---------------------------------------------------------------

if (length(out.lda) > 1){
  kvec <- sapply(out.lda, function(x) x@k)
  kchoose <- 50
  kchoose.i <- which(kvec == kchoose)
  out.lda <- out.lda[[kchoose.i]]
} else {
  kchoose <- out.lda@k
}

topics.mat <- posterior(out.lda)$topics
colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("CD41plus", cell))

m.umap.first <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = is.stem)) + geom_point(alpha = 0.5) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)

print(m.umap.first)

# Do variance calculation -------------------------------------------------

jfac <- 10^6
jpseudo <- 0
dat.impute.log <- log2(t(posterior(out.lda)$topics %*% posterior(out.lda)$terms) * jfac + jpseudo)
jchromos.num <- seq(25)
jchromos <- paste("chr", jchromos.num, sep = "")

cells.var.chromo.merged <- CalculateVarAll(dat.impute.log, jchromos)


# Plot with variance  -----------------------------------------------------

dat.umap.long.merge <- left_join(dat.umap.long, cells.var.chromo.merged)

m.umap.var <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_minimal() + 
  scale_colour_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_blank(), legend.position = "bottom") +
  xlab("") + ylab("") + 
  ggtitle(jmark)
print(m.umap.var)

m.umap.var.filt <- ggplot(dat.umap.long.merge %>% filter(umap2 < 2.5), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  scale_colour_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)
print(m.umap.var.filt)


# Get regions ------------------------------------------------------------

tm.result <- posterior(out.lda)

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


# Set up scRNA-seq WKM ----------------------------------------------------


# Load  -------------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.csv.gz"
inf.meta <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/tsne_clusterID_zebrafish_GateID_dataset.csv"
inf.meta2 <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/For_Jake/cluster_info_JY_edited.xlsx"

dat <- as.data.frame(fread(inf, stringsAsFactors = FALSE))
colnames(dat)[[1]] <- "gene"
rownames(dat) <- dat$gene
dat$gene <- NULL

meta <- fread(inf.meta, stringsAsFactors = TRUE)
meta2 <- read.xlsx2(inf.meta2, sheetIndex = 1, header = FALSE, colClasses = c("integer", "character")); colnames(meta2) <- c("ClusterID", "celltype")

colnames(meta) <- c("rname", "V1", "V2", "ClusterID", "experi")

meta <- left_join(meta, meta2)

# Get sum across celltypes  -----------------------------------------------

ctypes <- as.character(unique(meta$celltype))
names(ctypes) <- ctypes

ctype.sum <- lapply(ctypes, function(ctype){
  cells.tmp <- subset(meta, celltype == ctype)$rname
  cols.i <- which(colnames(dat) %in% cells.tmp)
  dat.tmp <- dat[, cols.i]
  genesums <- Matrix::rowSums(dat.tmp) 
  return(genesums)
})
genes <- names(ctype.sum[[1]])
assertthat::assert_that(all(sapply(ctype.sum, function(x) identical(names(x), genes))))

ctype.sum <- ctype.sum %>%
  bind_rows() %>%
  as.data.frame()

rownames(ctype.sum) <- genes

# boxplots
boxplot(log2(ctype.sum))

# normalize using DESeq2
ctype.metadata <- data.frame(celltype = colnames(ctype.sum))
rownames(ctype.metadata) <- colnames(ctype.sum)

dds <- DESeqDataSetFromMatrix(countData = round(ctype.sum), colData = ctype.metadata, design = ~1)
vsd <- vst(dds)

ctype.sum.vst <- assay(vsd)

pca.out <- prcomp(t(ctype.sum.vst), center = TRUE, scale. = TRUE, retx = TRUE)

plot(pca.out$x[, 1], pca.out$x[, 2], pch = 20)
text(pca.out$x[, 1], pca.out$x[, 2], labels = colnames(ctype.sum.vst))

boxplot(ctype.sum.vst)


# Remove lowly expressed genes? -------------------------------------------

xfilt <- 3.5
plot(density(unlist(ctype.sum.vst)))
abline(v = xfilt)

genes.keep.i <- which(apply(ctype.sum.vst, 1, function(x) quantile(x, 0.8)) > xfilt)
# genes.keep.i <- which(apply(ctype.sum.vst, 1, median) > xfilt)

ctype.sum.vst.filt <- ctype.sum.vst[genes.keep.i, ]

# do quantile normalization anyways
dat.mat.filt <- normalize.quantiles(ctype.sum.vst.filt, copy=TRUE)
colnames(dat.mat.filt) <- colnames(ctype.sum.vst.filt); rownames(dat.mat.filt) <- rownames(ctype.sum.vst.filt)

# dont do quantile normalization??
# dat.mat.filt <- ctype.sum.vst.filt

boxplot(dat.mat.filt)

jgenes <- sapply(rownames(dat.mat.filt), function(x) strsplit(x, "_")[[1]][[2]])

jgenes.dup <- jgenes[duplicated(jgenes)]
names(jgenes.dup) <- jgenes.dup

gene.exprs.collapsed <- lapply(jgenes.dup, function(jg){
  # get mean of duplicated genes
  apply(dat.mat.filt[which(jgenes == jg), ], 2, mean)
}) %>%
  do.call(rbind, .)

rows.exclude <- lapply(jgenes.dup, function(jg){
  which(jgenes == jg)
}) %>%
  unlist()

rows.include <- !seq(nrow(dat.mat.filt)) %in% rows.exclude
dat.mat.filt.collapsed <- dat.mat.filt[rows.include, ]

jgenes.collapsed <- sapply(rownames(dat.mat.filt.collapsed), function(x) strsplit(x, "_")[[1]][[2]])

# Get zscores -------------------------------------------------------------

dat.mat.filt.long <- tidyr::gather(data.frame(gene = jgenes.collapsed, dat.mat.filt.collapsed), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  ungroup()

dat.mat.filt.zscore <- tidyr::spread(dat.mat.filt.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>%
  as.data.frame()

rownames(dat.mat.filt.zscore) <- dat.mat.filt.zscore$gene

dat.mat.filt.zscore$gene <- NULL

# plot example genes
# jgene <- c("scpp8", "adam8a", "odc1", "lta4h", "thy1", "timp2b")
jgene <- c("rpl38", "med25", "pcdh7b", "pros1", "med25")
jgene <- c("mafba", "dapk1", "lhfpl6", "mrc1b")

ggplot(subset(dat.mat.filt.long, gene %in% jgene) %>% mutate(celltype = forcats::fct_reorder(celltype, dplyr::desc(zscore))), 
       aes(x = celltype, y = zscore)) +
  geom_point() + 
  geom_boxplot() + 
  ggtitle(jgene) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Look at topics and associate with celltypes -----------------------------

# order topics

topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
topics.sum$topic <- gsub("^X", "topic_", topics.sum$topic)

# show an island
topics.mat.tmp <- as.data.frame(topics.mat)
topics.mat.tmp$cell <- rownames(topics.mat)

topics.mat.tmp <- left_join(dat.umap.long, topics.mat.tmp)

keepns <- c(50, 100, 150, 200)
# keepn <- 25

for (keepn in keepns){
  print(paste("keepn", keepn))
  outpdf <- file.path(outdir, paste0("ZF_", jmark, "_topics_celltyping.keepn_", keepn, ".binarize_", jbin, ".kchoose_", kchoose, ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  print(m.umap.first)
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
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtop)
    
    print(m.umap)
    print(m.top)
    print(m.ctype)
  }
  dev.off()
}

# 
# # Project CD41 low cells --------------------------------------------------
# 
# inf.proj <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-ZFWKM-", 
#                    jmark, 
#                    "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/projections/", 
#                    "lda_out_meanfilt.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.OutObjs.K_30.x.ZF-ZFWKMCD41plus-", jmark, "_pcutoff_0.95_binfilt_cellfilt.2019-11-03.RData")
# assertthat::assert_that(file.exists(inf.proj))
# 
# load(inf.proj, v=T)
# 
# # do umap projected
# 
# dat.pred <- predict(umap.out, out.lda.predict$topics)
# 
# dat.pred.long <- data.frame(cell = rownames(dat.pred), umap1 = dat.pred[, 1], umap2 = dat.pred[, 2], stringsAsFactors = FALSE) %>%
#   mutate(is.stem = TRUE)
# 
# dat.umap.pred.merged <- bind_rows(dat.umap.long %>% mutate(is.stem = FALSE), dat.pred.long)
# 
# m.proj <- ggplot(dat.umap.pred.merged, aes(x = umap1, y = umap2, color = is.stem)) + geom_point(alpha = 0.5) +
#   theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_blank()) + 
#   xlab("") + ylab("")
# print(m.proj)
# 
# # Do louvain  -------------------------------------------------------------
# 
# dat.umap.long.louv <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long)
# 
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# m.louv <- ggplot(dat.umap.long.louv, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   ggtitle(jmark) + 
#   scale_color_manual(values = cbPalette)
# print(m.louv)
# 
# # plot outputs
# 
# pdf(file.path(outdir, paste0("ZF_downstream_", jmark, ".binarize_", jbin, ".kchoose_", kchoose, ".pdf")))
# 
#   print(m.umap.var)
#   print(m.umap.first)
#   print(m.louv)
#   print(m.proj)
# 
# dev.off()
# 
