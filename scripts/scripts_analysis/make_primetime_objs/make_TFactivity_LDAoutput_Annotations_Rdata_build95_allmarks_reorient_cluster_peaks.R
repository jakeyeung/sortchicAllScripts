# Jake Yeung
# Date of Creation: 2019-04-18
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_TFactivity_LDAoutput_Annotations_Rdata_build95_allmarks_reorient_cluster_peaks.R
# Analyze cluster peaks

rm(list=ls())

library(GGally)
library(purrr)

library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

library(igraph)

library(GGally)

library(biomaRt)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaData.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")

jpseudo <- 0
jscale <- 10^7

# Functions ---------------------------------------------------------------

PlotXYWithColor <- function(jsub, xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkblue"), jtitle = "", jcol.low = "gray85", jcol.mid = "gray50"){
  jrange <- range(jsub[[cname]])
  jmid <- min(jsub[[cname]]) + diff(range(jsub[[cname]])) / 2
  m1 <- ggplot(jsub, aes_string(x = xvar, y = yvar, color = cname)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_color_gradient2(low = "gray85", mid = "gray50", high = jcol, midpoint = jmid, limit = jrange) + 
    ggtitle(jtitle)
  return(m1)
}




# Load bulk sorted data ---------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Constants ---------------------------------------------------------------

jmark <- "H3K4me1"



# Load metadata -----------------------------------------------------------

experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
switch.rep.hash <- GetRepSwitchHash(experihash)

# Load UMAPs --------------------------------------------------------------

inf.spring <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"

load(inf.spring, v=T)

# Load TF activities ------------------------------------------------------

indir.mara <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_", 
                     jmark, 
                     ".filt_0.99.center_TRUE_K50")

mara.out <- LoadMARA(indir.mara, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = "")

# add activities to umap??

act.umap.long <- left_join(mara.out$act.long, dat.trajs.long %>% dplyr::filter(mark == jmark) %>% dplyr::select(cell, X1, X2))

# Plot top hits onto the UMAP ---------------------------------------------

jmotif <- "Foxc1"
jmotif <- "Pax6"

jmotif <- "Tal1"
jsub <- act.umap.long %>% filter(motif == jmotif) 
jmotif <- "Cebpb"

m1 <- PlotXYWithColor(jsub, jcol = scales::muted("darkblue"), jtitle = jmotif)
print(m1)
# jrange <- range(jsub$activity)
# jmid <- min(jsub$activity) + diff(range(jsub$activity)) / 2
# m1 <- ggplot(jsub, aes(x = X1, y = X2, color = activity)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_gradient2(low = "gray85", mid = "gray50", high = scales::muted("darkblue"), midpoint = jmid, limit = jrange) + 
#   ggtitle(jmotif)
# print(m1)




# Correlate with peak expression ------------------------------------------

# inf <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95.cells_from_bin_analysis/H3K4me1/mara_input/count_mats_peaks_norm/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50.txt"
# mat.exprs <- data.table::fread(inf, sep = "\t")

# load LDA
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95.cells_from_bin_analysis/", 
              jmark, 
              "/lda_input/lda_outputs.meanfilt_10.cellmin_0.cellmax_9999999.binarize.TRUE/lda_out_meanfilt.PZ-BM-", 
              jmark, ".CountThres0.K-25_50.Robj")
load(inf.lda, v=T)

out.lda <- out.lda[[length(out.lda)]]

docs.new <- SwitchColnames(unname(out.lda@documents))
out.lda@documents <- docs.new

tm.result <- posterior(out.lda)

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

# annotate only motifs
motifs.all <- unique(mara.out$zscores$motif)

mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')

gos <- getBM(motifs.all,
             attributes=c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand"),
             filters=c("external_gene_name"),
             mart=mart.obj) %>%
  filter(!startsWith(chromosome_name, prefix = "CHR")) %>%
  dplyr::mutate(chromosome_name = paste("chr", chromosome_name, sep = "")) %>%
  dplyr::rename(seqnames = chromosome_name,
                start = start_position,
                end = end_position) %>%
  mutate(strand = ifelse(strand == 1, "+", "-"))

gos.range <- makeGRangesFromDataFrame(gos)
gos.range$motif <- gos$external_gene_name

jnearest <- distanceToNearest(gos.range, regions.range)

motifs.query <- gos.range$motif[queryHits(jnearest)]
matches.query <- names(regions.range)[subjectHits(jnearest)]
distances.query <- mcols(jnearest)$distance

motif.peak.dist <- data.frame(motif = motifs.query, peak = matches.query, dist = distances.query, 
                              stringsAsFactors = FALSE)


# Plot motif and expression correlation -----------------------------------

mat.impute <- tm.result$topics %*% tm.result$terms 

mat.impute.filt <- mat.impute[, colnames(mat.impute) %in% motif.peak.dist$peak]

mat.impute.filt <- t(mat.impute.filt)

mat.impute.filt.long <- as.data.frame(mat.impute.filt) %>% mutate(peak = rownames(mat.impute.filt)) %>%
  tidyr::gather(key = "cell", value = "exprs", -peak)

mat.impute.filt.long <- left_join(mat.impute.filt.long, motif.peak.dist)

mat.impute.filt.long <- left_join(mat.impute.filt.long, act.umap.long)

mat.impute.filt.long$exprs <- log10(mat.impute.filt.long$exprs * jscale + jpseudo)

# plot top zscores

cors.all <- mat.impute.filt.long %>%
  group_by(motif) %>%
  summarise(pcor = cor(exprs, activity, method = "pearson")) %>%
  arrange(desc(abs(pcor)))

# add to zscores
zscores.cors <- left_join(mara.out$zscores, cors.all)

ggplot(cors.all, aes(x = pcor)) + geom_histogram() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle('correlations all')

# plot top 20 motifs
topn <- 20
# jmotif <- "Pitx2"
jmotifs <- (zscores.cors %>% arrange(desc(zscore)))$motif[1:topn]
jmotifs.bycor <- (zscores.cors %>% arrange(desc(abs(pcor))))$motif[1:topn]

pdf(paste0("~/data/scchic/pdfs/TF_analysis_gene_correlations.WithSpringTrajs.refmark.", jmark, ".", Sys.Date(), ".cluster_peaks_by_zscore.pdf"), useDingbats = FALSE)
for (jmotif in jmotifs){
  jsub <- mat.impute.filt.long %>% filter(motif == jmotif) 
  m1 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkred"), jtitle = jmotif)
  
  (jpeak <- subset(motif.peak.dist, motif == jmotif)$peak)
  (jdist <- subset(motif.peak.dist, motif == jmotif)$dist)
  (jcor <- signif(cor(jsub$activity, jsub$exprs), digits = 2))
  jzscore <- signif(subset(zscores.cors, motif == jmotif)$zscore, digits = 2)
  m2 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "exprs", jcol = scales::muted("darkblue"), jtitle = paste(jmotif, jpeak, jdist, jcor, jzscore))
  
  # plot correlations
  m.cor <- ggplot(subset(mat.impute.filt.long, motif == jmotif), aes(x = exprs, y = activity)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmotif, jpeak, jdist, jcor, jzscore)) + geom_smooth(method = "lm")
  multiplot(m1, m2, cols = 2)
  print(m.cor)
}
dev.off()

pdf(paste0("~/data/scchic/pdfs/TF_analysis_gene_correlations.WithSpringTrajs.refmark.", jmark, ".", Sys.Date(), ".cluster_peaks_by_cor.pdf"), useDingbats = FALSE)
for (jmotif in jmotifs.bycor){
  jsub <- mat.impute.filt.long %>% filter(motif == jmotif) 
  m1 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkred"), jtitle = jmotif)
  
  (jpeak <- subset(motif.peak.dist, motif == jmotif)$peak)
  (jdist <- subset(motif.peak.dist, motif == jmotif)$dist)
  (jcor <- signif(cor(jsub$activity, jsub$exprs), digits = 2))
  jzscore <- signif(subset(zscores.cors, motif == jmotif)$zscore, digits = 2)
  m2 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "exprs", jcol = scales::muted("darkblue"), jtitle = paste(jmotif, jpeak, jdist, jcor, jzscore))
  
  # plot correlations
  m.cor <- ggplot(subset(mat.impute.filt.long, motif == jmotif), aes(x = exprs, y = activity)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmotif, jpeak, jdist, jzscore)) + geom_smooth(method = "lm")
  
  multiplot(m1, m2, cols = 2)
  print(m.cor)
}
dev.off()


# Explore UMAP structure of the peak finding ------------------------------


custom.settings <- umap.defaults
custom.settings$n_neighbors <- 15
custom.settings$min_dist <- 0.4
custom.settings$metric <- "euclidean"

umap.out <- umap(tm.result$topics, config = custom.settings)

umap.out.long <- data.frame(cell = rownames(umap.out$layout),
                            umap1 = umap.out$layout[, 1],
                            umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)

betas.long <- as.data.frame(tm.result$topics) %>%
  mutate(cell = as.character(rownames(tm.result$topics))) %>%
  tidyr::gather(key = topic, value = beta, -cell)

ggplot(umap.out.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# explore tpics
# analyze topic matrix across cells
topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>%
  gather(key = "topic", value = "weight", -cell) %>%
  rowwise() %>%
  mutate(topic = as.numeric(substr(topic, 2, nchar(topic)))) %>%
  group_by(topic) %>%
  mutate(zscore = scale(weight, center = TRUE, scale = TRUE))

topics.sum <- topics.long %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(zscore < quantile(zscore, 0.97)) %>%
  mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)
print(topics.sum)

jtopics <- topics.sum$topic

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
barplot(height = topics.sum$entropy, names.arg = topics.sum$topic, xlab = "Topic", ylab = "Entropy Measure H", main = "Topics are now ordered by increasing entropy.", las = 2)


# add umap to this
umap.out.long <- left_join(betas.long, umap.out.long)

m.out <- ggplot(umap.out.long, aes(x = umap1, y = umap2, color = beta)) + geom_point() + facet_wrap(~topic) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.out)


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

# plot examples of peaks onto the UMAP?

top.peaks <- (top.peaks.annotated %>% filter(rnk < 2000))$term

# make long exprs

mat.impute.filt2 <- mat.impute[, colnames(mat.impute) %in% top.peaks]

mat.impute.filt2 <- t(mat.impute.filt2)

mat.impute.filt.long2 <- as.data.frame(mat.impute.filt2) %>% mutate(peak = rownames(mat.impute.filt2)) %>%
  tidyr::gather(key = "cell", value = "exprs", -peak) %>%
  mutate(exprs = log10(exprs * jscale + jpseudo))

# add UMAP information
mat.impute.filt.long2 <- left_join(mat.impute.filt.long2, dat.trajs.long %>% filter(mark == jmark) %>% dplyr::select(X1, X2, cell, louvain))


jmotif <- "Tal1"
# jgene2 <- "Cebpb"
jgene2 <- jmotif

# Compare with motifs
jsub <- mat.impute.filt.long %>% filter(motif == jmotif) 
m1 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "activity", jcol = scales::muted("darkred"), jtitle = jmotif)
(jpeak <- subset(motif.peak.dist, motif == jmotif)$peak)
(jdist <- subset(motif.peak.dist, peak == jpeak)$dist)
(jcor <- signif(cor(jsub$activity, jsub$exprs), digits = 2))
jzscore <- signif(subset(zscores.cors, motif == jmotif)$zscore, digits = 2)
m2 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "exprs", jcol = scales::muted("darkblue"), jtitle = paste(jmotif, jpeak, jdist, jcor, jzscore))
multiplot(m1, m2, cols = 2)



# (jpeak2 <- (top.peaks.annotated %>% ungroup() %>% filter(SYMBOL == jgene2) %>% filter(rnk == min(rnk)))$term)
jpeak2 <- jpeak
# take peaks
jsub2 <- mat.impute.filt.long2 %>% filter(peak == jpeak2)
print(head(jsub2))
m0 <- PlotXYWithColor(jsub2, xvar = "X1", yvar = "X2", cname = "exprs") + ggtitle(paste(jgene2, jpeak2))
print(m0)

multiplot(m0, m2, cols = 2)

# how sensitive are these peaks???
# scan along chromosome 2

jsub <- subset(mat.impute.filt.long2, grepl("^chr2:167[0-9]{6}", peak))

jsub$jstart <- as.numeric(sapply(jsub$peak, GetStart))
jsub$jend <- as.numeric(sapply(jsub$peak, GetEnd))
jsub$jchromo <- GetChromo(sapply(jsub$peak, GetChromo))

tss <- gsub(",", "", "chr2:167,688,913-167,689,428")
tss.start <- as.numeric(GetStart(tss))

jsub <- jsub %>% rowwise() %>% mutate(dist.to.cebpb = jstart - tss.start) %>% arrange(abs(dist.to.cebpb))
peaks.by.cebpb <- jsub[!duplicated(jsub$peak), ]

jpeaks <- peaks.by.cebpb$peak

pdf("~/data/scchic/pdfs/peaks_near_cebpb.pdf", useDingbats = FALSE)
for (jpeak2 in jpeaks){
  jgene2 <- subset(regions.annotated, region_coord == jpeak2)$SYMBOL
  jsub2 <- mat.impute.filt.long2 %>% filter(peak == jpeak2)
  jdist <- subset(peaks.by.cebpb, peak == jpeak2)$dist.to.cebpb
  m0 <- PlotXYWithColor(jsub2, xvar = "X1", yvar = "X2", cname = "exprs") + ggtitle(paste(jgene2, jpeak2, jdist))
  print(m0)
}
dev.off()



#
# 
# # Layer on bulk expression data -------------------------------------------
# 
# pdf(paste0("~/data/scchic/pdfs/lda_cluster_peak_output.", Sys.Date(), ".pdf"), useDingbats = FALSE)
# topn <- 150
# topn.plot <- 50
# for (jtopic in jtopics){
#   top.genes <- subset(top.peaks.annotated, topic == jtopic)$SYMBOL[1:topn]
#   jsub <- subset(dat.long, Gene_Name %in% top.genes)
#   jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
#   jlevels <- as.character(jsub.sorted.summarised$CellType)
#   jsub$CellType <- factor(jsub$CellType, levels = jlevels)
#   
#   
#   m.out <- ggplot(umap.out.long %>% filter(topic == jtopic), aes(x = umap1, y = umap2, color = beta)) + geom_point() + facet_wrap(~topic) +
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#   print(m.out)
#   
#   m.ctype <- ggplot(jsub,
#                     aes(x = CellType , y = zscore)) +
#     geom_boxplot() +
#     # geom_violin() +
#     geom_jitter(width = 0.1, size = 0.5) +
#     # geom_line() +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle(jtopic)
#   print(m.ctype)
#   
#   
#   
#   # plot top bins?
#   top.sub <- top.peaks.annotated %>% filter(topic == jtopic & rnk < topn.plot) %>% arrange(rnk)
#   top.sub <- OrderDecreasing(top.sub, jfactor = "term", jval = "beta")
#   m.tops <- ggplot(top.sub, aes(x = term, y = beta, label = SYMBOL)) +
#     geom_point() +
#     geom_text_repel() +
#     theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(),
#                        panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ggtitle(paste("Topic:", jtopic))
#   print(m.tops)
# }
# 
# dev.off()
# 
# 
