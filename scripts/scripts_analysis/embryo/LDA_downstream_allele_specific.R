# Jake Yeung
# Date of Creation: 2019-03-28
# File: ~/projects/scchic/scripts/scripts_analysis/embryo/LDA_downstream_allele_specific.R
# Allele-specific outputs??

rm(list=ls())

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


source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")
source("scripts/Rfunctions/PlotFunctions.R")


SwapChromoFromCoord <- function(coord.old, chromo.hash, jsep = ":", outsep = "_"){
  # chr1:3000000-3100000 -> chr1_b6:2000000-31000000
  jsplit <- strsplit(coord.old, split = jsep)[[1]]
  chromo.old <- jsplit[[1]]
  coord.startend <- paste(jsplit[2:length(jsplit)], collapse = jsep)
  chromo.new <- chromo.hash[[chromo.old]]
  assertthat::assert_that(!is.null(chromo.new))
  coord.new <- paste(chromo.new, coord.startend, sep = jsep)
  return(coord.new)
}

StripSuffixFromCoord <- function(coord.full){
  # chrX_b6:48500000-48600000 -> chrX:48500000-48600000
  jsplit <- strsplit(coord.full, ":")[[1]]
  chromo <- jsplit[[1]]
  startend <- jsplit[[2]]
  chromo.new <- strsplit(chromo, "_")[[1]][[1]]
  return(paste0(chromo.new, ":", startend))
}

GetChromoSplitAllele <- function(x){
  xtmp <- GetChromo(x)
  return(strsplit(xtmp, "_")[[1]][[1]])
}

# Load stuff --------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/from_cluster/embryo/lda_out_meanfilt.mouse_embryo_K36me3_build95_AS_LDA25.Robj"

load(inf, v=T)

out.lda <- out.lda[[1]]
kchoose <- out.lda@k



# Wrangle terms to Bcl6 vs Castaneous -------------------------------------

nchromos <- 20
chrs.raw <- seq(2 * nchromos)
chrs.new <- chrs %% nchromos  # chromosome 0 is X chromosome
chrs.allele <- ceiling(chrs.raw / nchromos)
chrs.allele.named <- sapply(chrs.allele, function(x) ifelse(x == 1, "b6", "cast"))

chrs.dat <- data.frame(chromo.raw = paste("chr", chrs.raw, sep = ""),
                       chromo.num = as.character(chrs.raw), 
                       chromo.mod = paste("chr", chrs.new, sep = ""),
                       allele = chrs.allele.named, 
                       stringsAsFactors = FALSE)
chrs.dat$chromo.mod <- sapply(chrs.dat$chromo.mod, function(x) ifelse(x == "chr0", "chrX", x))
chrs.dat$chromo.full <- paste(chrs.dat$chromo.mod, chrs.dat$allele, sep = "_")

chromo.mod.hash <- hash::hash(chrs.dat$chromo.raw, chrs.dat$chromo.mod)
chromo.full.hash <- hash::hash(chrs.dat$chromo.raw, chrs.dat$chromo.full)
chromo.num.hash <- hash::hash(chrs.dat$chromo.num, chrs.dat$chromo.full)

# rename terms
terms.new <- sapply(out.lda@terms, SwapChromoFromCoord, chromo.full.hash, jsep = ":", outsep = "_")
terms.names.new <- sapply(names(out.lda@terms), SwapChromoFromCoord, chromo.num.hash, jsep = "_", outsep = "-")
names(terms.new) <- terms.names.new

# swap values of out.lda@terms
out.lda@terms <- terms.new


# Process data ------------------------------------------------------------

tm.result <- posterior(out.lda)


nn=30
jmetric='euclidean' 
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)

dat.umap <- umap(tm.result$topics, config = custom.settings)


# Plot output -------------------------------------------------------------

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))

m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, tm.result$topics, c("pink", "darkred"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)

plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, asp = 0.75, cex = 0.2)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1), pty = "s")
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75, cex = 0.2)
}, jcol.rgbs, seq(kchoose))


# Sort topics by interestingness ------------------------------------------

# analyze topic matrix across cells
topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>% 
  gather(key = "topic", value = "weight", -cell) %>%
  rowwise() %>%
  mutate(topic = as.numeric(substr(topic, 2, nchar(topic)))) %>%
  group_by(topic) %>%
  mutate(zscore = scale(weight, center = TRUE, scale = TRUE))

topics.sum <- topics.long %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(zscore < quantile(zscore, 0.99)) %>%
  mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)
print(topics.sum)

barplot(height = topics.sum$entropy, names.arg = topics.sum$topic, xlab = "Topic", ylab = "Entropy Measure H", main = "Low H suggests topic with islands of high cell weights. High H evenly distributed weights")

ggplot(topics.long, aes(x = topic, y = zscore, group = topic)) + geom_boxplot() + theme_bw() + scale_x_continuous(breaks = seq(out.lda@k))

# plot distributions
ggplot(topics.long, aes(x = zscore, group = topic)) + geom_density() + facet_wrap(~topic) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(topics.long, aes(x = zscore, group = topic)) + geom_histogram(bins = 50) + facet_wrap(~topic) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(topics.long, aes(x = zscore, group = topic)) + geom_histogram(bins = 50) + facet_wrap(~topic) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_log10()

# summarize with a number???
# ggplot(topics.sum, aes(x = topic, y = entropy)) + geom_bar(stat = "identity") + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggplot(topics.sum, aes(x = entropy)) + geom_density() 


# Plot interesting topics based on entropy  -------------------------------

topics.ranked <- topics.sum$topic


# Summarize terms  --------------------------------------------------------


top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))
top.peaks$coord <- sapply(top.peaks$term, StripSuffixFromCoord)
regions <- data.frame(seqnames = sapply(colnames(tm.result$term), GetChromoSplitAllele),
                      start = sapply(colnames(tm.result$term), GetStart),
                      end = sapply(colnames(tm.result$term), GetEnd), 
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$term)

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)
top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))
