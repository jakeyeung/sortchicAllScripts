# Jake Yeung
# Date of Creation: 2019-02-19
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/figures_for_presentation.R
# Figure for lab meeting


rm(list=ls())

jstart <- Sys.time()

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(umap)
library(ggrepel)
library(biomaRt)
library(igraph)  # louvain
library(Gviz)
library(GenomicRanges)

library(Cairo)
library(ggrastr)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/MaraDownstream.R")

# Load LDA data -----------------------------------------------------------


jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
jcolvec <- c("blue", "gray80", "red")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
names(out.objs) <- jmarks
names(out.objs.nobin) <- jmarks

tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))

# get count matrices
count.mat.lst <- lapply(out.objs.nobin, function(x) x$count.mat)

# what's the distirbution of counts?
lapply(count.mat.lst, range)

# sparsity of matrix?

sparsity <- lapply(count.mat.lst, function(count.mat) 1 - nnzero(count.mat) / length(count.mat))
multicounts <- lapply(count.mat.lst, function(count.mat) length(which(as.matrix(count.mat) > 1))  / length(count.mat))

plot(density(unlist(as.matrix(count.mat.lst[[3]]))))

# cells with large counts, what are they lonely ?
# indx <- which(count.mat.lst[[1]] > 50, arr.ind = TRUE)

# Visualize sparsity of matrix --------------------------------------------

jseed <- 0
jsize <- 0.01
i <- 1
jsub <- count.mat.lst[[i]][sample(x = seq(nrow(count.mat.lst[[i]])), size = jsize * nrow(count.mat.lst[[i]]), replace = FALSE), ]
# order cells by decreasing size?
jsub <- jsub[order(Matrix::rowSums(jsub), decreasing = TRUE), order(Matrix::colSums(jsub), decreasing=TRUE)]

pdf("~/Documents/presentations_postdoc/2019-02-20_meeting.pdf", useDingbats = FALSE)

image(jsub, xlab = "Cells", ylab = "Bins", main = paste(jmarks[[i]]), useRaster = TRUE)


# Summarize mean and variance relationships -------------------------------

mean.vars.lst <- lapply(jmarks, function(jmark){
  cell.means <- Matrix::colMeans(count.mat.lst[[jmark]])
  cell.vars <- apply(count.mat.lst[[jmark]], 2, var)
  cell.cv2 <- cell.vars / cell.means ^ 2
  dat.tmp <- data.frame(mark = jmark, cell.means = cell.means, cell.cv2 = cell.cv2)
  return(dat.tmp) 
})

mean.vars.dat <- bind_rows(mean.vars.lst)

m.meanvar <- ggplot(mean.vars.dat, aes(x = log10(cell.means),  y = log10(cell.cv2))) +  
  theme_bw(20) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point_rast(size = 1, alpha = 0.25) + facet_wrap(~mark) + 
  geom_abline(slope = -1, intercept = 0) + 
  xlab("log10(Cell Mean)") + ylab("log10(CV2)")
print(m.meanvar)

dev.off()


# Show the cell types for H3K4me1 -----------------------------------------

jmark <- "H3K4me1"
tm.result <- posterior(out.objs[[jmark]]$out.lda)
topics.mat <- tm.result$topics
topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)

jcolvec <- c(muted("blue"), "gray80", muted("red"))
nn <- 40
jmetric <- "euclidean"
jmindist <- 0.2
jseed <- 123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)

dat.umap <- umap(topics.mat, config = custom.settings)

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)

jtopic <- 1
m <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = paste0("`", jtopic, "`"))) + geom_point() +
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_continuous(paste0("Topic ", jtopic, "\nWeight")) +
  scale_color_gradient2(name = paste0("Topic ", jtopic, "\nWeight"),
                        low = jcolvec[[1]], mid = jcolvec[[2]], high = jcolvec[[3]], 
                        midpoint = quantile(range(dat.umap.long[[as.character(jtopic)]]), 0.3))
print(m)

# show top hits for topic
top.peaks <- tidytext::tidy(out.objs[[jmark]]$out.lda, matrix = "beta") %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd),
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)
top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))

i <- jtopic
m.top <- subset(top.peaks.annotated, topic == i) %>% top_n(n=25, wt = beta) %>%
  mutate(term = forcats::fct_reorder(term, dplyr::desc(beta))) %>%
  ggplot(aes(x = term, y = log10(beta), label = SYMBOL)) +
  geom_point() + theme_bw(20) + geom_text_repel(size = 8) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Log10 Topic-to-Region Weight") +
  ggtitle(paste("Top peak weights for topic:", i))
print(m.top)

# show cell-type specific expression on chart on these 25

dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

topn <- 25
top.genes <- subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]
# top.genes <- "Lancl2"
jsub <- subset(dat.long, Gene_Name %in% top.genes)
jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
jlevels <- as.character(jsub.sorted.summarised$CellType)
jsub$CellType <- factor(jsub$CellType, levels = jlevels)
m.celltype <- ggplot(jsub, 
       aes(x = CellType , y = zscore)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_classic(20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(paste("Genes assigned to top 25 peaks in topic", jtopic)) + xlab("")
print(m.celltype)

# show Hox genes
jsub <- subset(dat.long, grepl("Hox", Gene_Name))
jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
jlevels <- as.character(jsub.sorted.summarised$CellType)
jsub$CellType <- factor(jsub$CellType, levels = jlevels)
m.celltype.hox <- ggplot(jsub, 
                     aes(x = CellType , y = zscore)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_classic(15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle(paste("Expression of Hox genes")) + xlab("")
print(m.celltype.hox)

# Louvain  ----------------------------------------------------------------

custom.settings2 <- custom.settings
custom.settings2$n_neighbors <- 27
dat.umap.forlouvain <- umap(topics.mat, config = custom.settings2)
cell.indx <- hash(rownames(dat.umap.forlouvain$knn$indexes), dat.umap.forlouvain$knn$indexes[, 1])
cell.indx.rev <- hash(dat.umap.forlouvain$knn$indexes[, 1], rownames(dat.umap.forlouvain$knn$indexes))
nr <- nrow(dat.umap.forlouvain$knn$indexes)
nc <- ncol(dat.umap.forlouvain$knn$indexes)
edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
colnames(edgelist) <- c("from", "to")
for (vertex.i in seq(nr)){
  istart <- nc*(vertex.i - 1)+1
  iend <- nc*vertex.i
  edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
  edgelist[istart : iend, 2] <- sapply(dat.umap.forlouvain$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
  # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
}
g <- graph_from_data_frame(edgelist, directed=FALSE)
g.out <- cluster_louvain(g, weights = NULL)
V(g)$color <- g.out$membership
clstr <- hash(g.out$names, g.out$membership)

dat.umap.long$louvain <- sapply(dat.umap.long$cell, function(x) clstr[[x]])

m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) + 
  geom_point(size = 0.8) +
  theme_bw(12) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "bottom") + guides(color=guide_legend(title="Cluster"))
print(m.louvain)



# Look at H3K27me3 --------------------------------------------------------

jmark <- "H3K27me3"
tm.result <- posterior(out.objs[[jmark]]$out.lda)
topics.mat <- tm.result$topics
topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)

jcolvec <- c(muted("blue"), "gray80", muted("red"))
nn <- 40
jmetric <- "euclidean"
jmindist <- 0.4
jseed <- 123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)

dat.umap <- umap(topics.mat, config = custom.settings)

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)


# show top hits for topic
top.peaks <- tidytext::tidy(out.objs[[jmark]]$out.lda, matrix = "beta") %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd),
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)
top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))


topn <- 500
pdf("~/Documents/presentations_postdoc/H3K27me3_topic_colors.pdf", useDingbats = FALSE)
jtopics <- seq(ncol(topics.mat))
for (jtopic in jtopics){
  m <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = paste0("`", jtopic, "`"))) + geom_point() +
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # scale_color_continuous(paste0("Topic ", jtopic, "\nWeight")) +
    scale_color_gradient2(name = paste0("Topic ", jtopic, "\nWeight"),
                          low = jcolvec[[1]], mid = jcolvec[[2]], high = jcolvec[[3]], 
                          midpoint = quantile(range(dat.umap.long[[as.character(jtopic)]]), 0.3))
  m.top <- subset(top.peaks.annotated, topic == jtopic) %>% top_n(n=50, wt = beta) %>%
    mutate(term = forcats::fct_reorder(term, dplyr::desc(beta))) %>%
    ggplot(aes(x = term, y = log10(beta), label = SYMBOL)) +
    geom_point() + theme_bw(12) + geom_text_repel(size = 5) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("Log10 Topic-to-Region Weight") +
    ggtitle(paste("Top peak weights for topic:", jtopic))
  
  top.genes <- subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]
  jsub <- subset(dat.long, Gene_Name %in% top.genes)
  jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
  jlevels <- as.character(jsub.sorted.summarised$CellType)
  jsub$CellType <- factor(jsub$CellType, levels = jlevels)
  m.celltype <- ggplot(jsub, 
                       aes(x = CellType , y = zscore)) + 
    geom_boxplot() +
    geom_jitter(width = 0.1, size = 0.5) +
    theme_classic(20) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(paste("Top", topn, "peaks in topic", jtopic)) + xlab("")
  print(m.celltype)
  print(m)
  print(m.top)
}
dev.off()


# Do motif analysis -------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3")

mdirs <- lapply(jmarks, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.repress <- c("H3K27me3", "H3K9me3")
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  print(mdir)
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all

mdirs <- c(mdirs, mdirs.repress)
mara.outs <- lapply(mdirs, LoadMARA)

act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long, mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)


zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, mara.outs[[3]]$zscores, mara.outs[[4]]$zscores), left_join, by = "motif")
cnames <- c("motif", paste("zscore", jmarks.all, sep = "."))
colnames(zscores.merged) <- cnames
zscore.thres <- 0.75
zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
  ifelse(max(jrow[[2]], jrow[[3]]) > zscore.thres, jrow[[1]], NA)
})
zscores.merged.mean <- zscores.merged
zscores.merged.mean$zscore.mean <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, mean)
zscores.merged.mean$zscore.max <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, max)

zscores.merged.mean <- zscores.merged.mean %>%
  arrange(desc(zscore.mean))

act.mat.merged <- spread(act.long.merged, key = "cell", value = "activity")


# UMAP on the terms
motifs.filt2 <- subset(zscores.merged.mean, !is.na(motif.lab))$motif.lab
motifs.filt <- zscores.merged$motif
jsub <- subset(act.mat.merged, motif %in% motifs.filt)
# umap.out <- umap(t(subset(jsub, select = -motif)))  # maybe interesting? need to color by batch

nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.4, 0.2, 0.1)
jmetric <- "euclidean"
jseed=123
custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))

out.lda.lst <- lapply(jmarks.all, function(jmark) out.objs[[jmark]]$out.lda)

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

names(umap.lda.lst) <- jmarks.all
names(mara.outs) <- jmarks.all

dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

jmotif <- "Tal1"

jmotifs <- c("Cebpb", "Tal1", "Pml", "Nfatc1", "Ebf1", "Bcl3", "Bptf", "Stat6")

pdf("~/Documents/presentations_postdoc/motif_analysis.pdf", useDingbats = FALSE)
for (jmotif in jmotifs){
  plts.lst <- lapply(jmarks.all, function(jmark){
    return(PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]], mara.outs[[jmark]]$zscores, jmark, jsize = 0.5))
  })
  multiplot(plts.lst[[1]], plts.lst[[2]], cols = 2)
  multiplot(plts.lst[[1]], plts.lst[[3]], plts.lst[[2]], plts.lst[[4]], cols = 2)
  multiplot(plts.lst[[3]], plts.lst[[4]], cols = 2)
  
}
dev.off()
print(Sys.time() - jstart)

