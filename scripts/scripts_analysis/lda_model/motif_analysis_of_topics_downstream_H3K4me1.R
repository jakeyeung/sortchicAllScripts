# Jake Yeung
# Date of Creation: 2019-01-21
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/motif_analysis_of_topics_downstream.R
# 
# 

rm(list=ls())

library(dplyr)
library(topicmodels)
library(RcisTarget)
library(feather)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(JFuncs)
library(umap)
library(ggplot2)
library(forcats)
library(ggrepel)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

# Functions ---------------------------------------------------------------

LoadLDA <- function(jchip){
  jdist <- 1000L
  jmean <- 1
  jmin <- 100L
  jmax <- 500000L
  # binarize <- "TRUE";  jtops <- "15_20_25_30_35"
  binarize <- "FALSE"; jtops <- "15_20_25_30_35"
  
  jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
  inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
  infbase <- basename(inf)
  infbase <- strsplit(infbase, ".Robj")[[1]][[1]]
  
  # 0.98 threshold 
  # inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
  # 0.96 threshold 
  inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.0.96.Robj"))
  
  load(inf, v=T)
  load(inf.GREAT, v=T)
  
  if (exists("regions.range")){
    regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                    TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                    annoDb='org.Mm.eg.db'))
    regions.annotated$region_coord <- rownames(regions.range)
  }
  tm.result <- posterior(out.lda)
  return(list('out.lda' = out.lda, 'tm.result' = tm.result, 'topic.regions' = topic.regions, 'regions.annotated' = regions.annotated))
}

PlotImputedGenes <- function(tm.result, regions.annotated, jgene, jdist, jchip, show.plot=TRUE, return.plot.only=FALSE){
  regions.sub <- subset(regions.annotated, abs(distanceToTSS) < jdist & grepl(jgene, SYMBOL)) %>%
    mutate(coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"),
           coord.kb = paste(paste(seqnames, start / 1000, sep = ":"), end / 1000, sep = "-"))
  if (nrow(regions.sub) == 0){
    warning(paste("empty regions.sub for jgene", jgene, "jchip", jchip))
    m <- ggplot()
    return(m)
  }
  peaks.keep <- regions.sub$coord
  peaks.keep.kb <- regions.sub$coord.kb
  
  mat.norm <- t(tm.result$topics %*% tm.result$terms)
  
  jlab <- jgene
  
  # run on imputed matrix
  row.i <- which(rownames(mat.norm) %in% peaks.keep)
  
  if (length(row.i) > 1){
    # print(paste("Merging", length(row.i)), "rows")
    jcounts.norm <- Matrix::colSums(mat.norm[row.i, ])
  } else {
    jcounts.norm <- mat.norm[row.i, ]
  }
  # print(range(jcounts.norm))
  jcol.counts <- ColorsByCounts(jcounts.norm, nbreaks = 100, colvec = c("lightblue", "darkblue"))
  # plot topics soft clustering weights
  topics.mat <- tm.result$topics
  nn=5
  nnterms=15
  jmetric='euclidean' 
  jmindist=0.1
  custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
  custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
  dat.umap <- umap(topics.mat, config = custom.settings)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  jpeaks.str <- paste(peaks.keep.kb, collapse = ",")
  jmain <- paste0(jchip, " ", jlab, ' npeaks ', length(row.i), "\n", jpeaks.str)
  
  # prepare plot object
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2], 
                    jcol = jcol.counts,
                    counts.norm = jcounts.norm)
  m <- ggplot(dat, aes(x = umap1, y = umap2, col = jcol.counts)) + geom_point() + scale_color_identity() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    ggtitle(jmain)
  
  if (show.plot){
    print(m)
  }
  
  # plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, 
  #      main = jmain, col = jcol.counts, asp = 0.75)
  if (!return.plot.only){
    return(list('regions.sub' = regions.sub, 'm' = m))
  } else {
    return(m)
  }
}

PlotAllMarks <- function(jgene, jmarks, jdist, out.objs){
  out <- lapply(jmarks, function(jchip){
    regions.annotated <- out.objs[[jchip]]$regions.annotated
    tm.result <- out.objs[[jchip]]$tm.result
    m <- PlotImputedGenes(tm.result, regions.annotated, jgene, jdist, jchip, show.plot=FALSE, return.plot.only = TRUE)
  })
  return(out)
}

# Load a topic output -----------------------------------------------------

# jchip <- "H3K4me1"
jchips <- c("H3K27me3", "H3K4me1", "H3K4me3", "H3K9me3")

jchip <- jchips[[2]]

etables <- readRDS(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/motif_enrichment_tables.", jchip, ".rds"))

# load the LDA output for visualization

jchips <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
out.objs <- lapply(jchips, LoadLDA)
names(out.objs) <- jchips

etables[[17]][1:15, 1:5]  # Tal1
etables[[12]][1:15, 1:5]
etables[[19]][1:15, 1:5]

etables.all <- bind_rows(lapply(seq(length(etables)), function(i) etables[[i]] %>% mutate(topic = i))) %>%
  arrange(desc(NES))

# Show the topic plots to remind  -----------------------------------------

out.lda <- out.objs[[jchip]]$out.lda

(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)
topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

nn=40
nnterms <- 15
jmetric='euclidean' 
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)

topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)
dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)


jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))


# Annotation of peaks -----------------------------------------------------

plot(density(out.objs[["H3K4me1"]]$regions.annotated$width))
plot(density(out.objs[["H3K4me3"]]$regions.annotated$width))

plot(density(log10(abs(out.objs[["H3K4me1"]]$regions.annotated$distanceToTSS) + 1)))
plot(density(log10(abs(out.objs[["H3K4me3"]]$regions.annotated$distanceToTSS) + 1)))

regions.annot.merged <- bind_rows(
  out.objs[["H3K4me1"]]$regions.annotated %>% mutate(mark = "H3K4me1"),
  out.objs[["H3K4me3"]]$regions.annotated %>% mutate(mark = "H3K4me3")
)

ggplot(regions.annot.merged, aes(x = mark, y = log10(abs(distanceToTSS) + 1))) + 
  geom_violin() + geom_boxplot() + theme_bw()
# how many at promoter?

regions.annot.merged %>%
  rowwise() %>%
  mutate(annotshort = strsplit(annotation, " ")[[1]][[1]]) %>%
  group_by(annotshort, mark) %>%
  summarise(ann.count = length(annotation)) %>%
  group_by(mark) %>%
  mutate(ann.frac = ann.count / sum(ann.count)) %>%
  ggplot(aes(x = annotshort, y = ann.frac, group = mark, fill = mark)) + 
  geom_bar(stat = "identity", position = position_dodge()) + theme_bw() +
  ggtitle("H3K4me1 Peak Annotations") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylim(c(0, 0.5)) + 
  xlab("") + ylab("Fraction of annotation")


# Get top hits ------------------------------------------------------------

top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))
# annotate regions?
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



# Motif analysis ----------------------------------------------------------

# show topic 17: Hbb and Sox6 enriched here 
jtopic <- 17
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `17`)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
print(m)

# which motifs here?
peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
out.objs[["H3K4me1"]]$regions.annotated

print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 4 & AUC > 0.005))



# Do for all topics -------------------------------------------------------

outdir <- paste0("~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/", jchip)
dir.create(outdir)
pdf(file.path(outdir, paste0("topic_weights_from_peak_calling.", jchip, ".pdf")), useDingbats = FALSE)


par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))

# show top hits of betas
for (i in seq(kchoose)){
  # etables at topic 17
  m.topmotifs <- etables[[i]] %>% top_n(n = 25, wt = NES) %>% 
    mutate(motif = forcats::fct_reorder(motif, dplyr::desc(NES))) %>%
    rowwise() %>%
    mutate(TFname = sapply(TF_highConf, function(x) strsplit(x, " ")[[1]][1])) %>%
    ggplot(aes(x = motif, y = NES, label = TFname)) + 
    geom_text_repel() +
    geom_point() + 
    theme_bw(12) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab("") + ggtitle(paste("Motif enrichment topic:", i))
  print(m.topmotifs)
}



jtopics <- seq(kchoose)
for (jtopic in jtopics){
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  # handle variable column names
  # https://stackoverflow.com/questions/13445435/ggplot2-aes-string-fails-to-handle-names-starting-with-numbers-or-containing-s
  m <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2", color = paste("`", jtopic, "`", sep=""))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
  print(m)
  # which motifs here?
  peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
  regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
  etab.out <- subset(etables[[jtopic]], select = -enrichedGenes, NES > 3.5 & AUC > 0.005)
  write.table(etab.out, file.path(outdir, paste0(jchip, ".topic_", jtopic, "_topmotifs.txt")), quote = FALSE, sep = "\t", row.names = FALSE)
}
dev.off()

# 
# # Show topic 1 ------------------------------------------------------------
# 
# jtopic <- 1
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `1`)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
# print(m)
# 
# # which motifs here?
# peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
# regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
# # out.objs[["H3K4me1"]]$regions.annotated
# 
# # print(etables[[jtopic]]$enrichedGenes[[1]])
# # print(subset(etables[[jtopic]], select = -enrichedGenes)[1:10, ])
# print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 4 & AUC > 0.005))
# 
# # Show topic 12 ------------------------------------------------------------
# 
# jtopic <- 12
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `12`)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
# print(m)
# 
# # which motifs here?
# peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
# regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
# 
# print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 4 & AUC > 0.005))
# 
# 
# # Show tpic 19 ------------------------------------------------------------
# 
# jtopic <- 19
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `19`)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
# print(m)
# 
# # which motifs here?
# peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
# regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
# 
# print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 3.5 & AUC > 0.005))
# 
# # Show tpic 10 ------------------------------------------------------------
# 
# jtopic <- 10
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `10`)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
# print(m)
# 
# # which motifs here?
# peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
# regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
# 
# print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 3.5 & AUC > 0.005))
# 
# # Show tpic 9 ------------------------------------------------------------
# 
# jtopic <- 9
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `9`)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
# print(m)
# 
# # which motifs here?
# peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
# regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
# 
# print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 3.5 & AUC > 0.005))
# 
# # Show tpic 20 ------------------------------------------------------------
# 
# jtopic <- 20
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = `9`)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
# print(m)
# 
# # which motifs here?
# peaks.sub <- out.objs[[jchip]]$topic.regions[[jtopic]]
# regions.sub <- subset(regions.annotated, abs(distanceToTSS) < 10000 & region_coord %in% peaks.sub)
# 
# print(subset(etables[[jtopic]], select = -enrichedGenes, NES > 3.5 & AUC > 0.005))
# 
# 
# 
# 
# # Tables ------------------------------------------------------------------
# 
# # look at topic 1 H3K4me1
# # Spib motif is upregulated
# # Target genes: Cebpg
# 
# print(etables[[1]]$enrichedGenes[[1]])
# jgene <- "Cepbg"
# jgene <- "Sox6"
# jgene <- "Igf2r"
# jdist <- 10000
# system.time(
#   plots.out <- PlotAllMarks(jgene, jchips, jdist, out.objs)
# )
# multiplot(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], cols = 4)
# 
# # Sox10 outputs
# 
# head(etables[[3]]$TF_highConf, n = 30)
# etables[[3]]$enrichedGenes[[1]]
# 
# jgene <- "Dgkg"
# jdist <- 10000
# system.time(
#   plots.out <- PlotAllMarks(jgene, jchips, jdist, out.objs)
# )
# multiplot(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], cols = 4)
# 
# # Irf?
# head(etables[[3]]$TF_highConf, n = 30)
# etables[[3]]$enrichedGenes[[6]]
# 
# jgene <- "Irf4"
# system.time(
#   plots.out <- PlotAllMarks(jgene, jchips, jdist, out.objs)
# )
# multiplot(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], cols = 4)
