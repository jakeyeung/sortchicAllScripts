# Jake Yeung
# Date of Creation: 2019-02-04
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/plot_marker_genes_all_4_marks.R
# Plot marker genes on all 4 markers on bins analysis

library(dplyr)
library(ggplot2)
library(JFuncs)
library(topicmodels)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(hash)
library(umap)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")


# load from bins ----------------------------------------------------------

jsize <- 0.5
jcolvec <- c("blue", "yellow", "red")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

out.objs <- lapply(jmarks, LoadLDABins)
names(out.objs) <- jmarks

tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))

# Plot hit across genes ---------------------------------------------------

# settings for UMAP
# single setting for all 4 
nn=40
# nn=15
nnterms <- 15
jmetric='euclidean'
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)

regions.annot <- out.objs[[1]]$regions.annotated

jgene <- "Sox6"


jgene <- "Tal1"
jsub <- subset(regions.annot, grepl(jgene, SYMBOL)) %>% arrange(abs(distanceToTSS))


# pick peak with largest weight for HK4me1
jpeaks <- jsub$region_coord
jmark <- "H3K4me1"
# jmark <- "H3K27me3"

jmax <- apply(tm.result.lst[[jmark]]$terms[, jpeaks], 2, max)
jmax.i <- apply(tm.result.lst[[jmark]]$terms[, jpeaks], 2, which.max)  # which topic max occurs. Often they agree.
jpeak <- jpeaks[which(jmax == max(jmax))]

which.max(tm.result.lst[[jmark]]$terms[, jpeaks])

# find neutrophilmarkers
ref.mark <- "H3K4me1"
jgene <- "S100a8"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)

ref.mark <- "H3K27me3"
jgene <- "Hoxc13"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
(jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]]))
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], 
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
m13 <- m.lst[[3]]

ref.mark <- "H3K27me3"
jgene <- "Hoxc4"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
(jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]]))
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], 
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
m4 <- m.lst[[3]]

jpeak13 <- "chr15:102820000-102920000"
# jpeak4 <- "chr15:103020000-103120000"
jpeak4 <- "chr5:139540000-139640000"  # copy and pasting dashes can be a bug!!
# jpeak4.works <- "chr5:139540000-139640000"
mat.norm <- t(tm.result.lst[["H3K27me3"]]$topics %*% tm.result.lst[["H3K27me3"]]$terms)
mat.sub <- mat.norm[c(jpeak13, jpeak4), ]
# mat.sub <- mat.norm[c(jpeak13, jpeak4.works), ]
mat.sub.max <- apply(mat.sub, 1, max)

# they look identical??
plot(mat.sub[1, ], mat.sub[2, ])
plot(mat.sub[1, ], mat.sub[2, ], log = "xy")

terms.sub <- t(tm.result.lst[["H3K27me3"]]$terms[, c(jpeak13, jpeak4)])
plot(log(terms.sub[1, ]), log(terms.sub[2, ]))

multiplot(m13, m4)

# plot whole Hoxc cluster, 4 to 13
ref.mark <- "H3K27me3"
jgrep <- "Hoxc"
jsub <- subset(out.objs[[ref.mark]]$regions.annot, grepl("Hoxc", SYMBOL)) %>% 
  rowwise() %>%
  mutate(suffix = as.numeric(strsplit(SYMBOL, jgrep)[[1]][[2]])) %>%  # get number from string
  arrange(suffix)
jgenes <- unique(jsub$SYMBOL)
for (jgene in jgenes){
  out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
  jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]])
  m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], 
                                                           show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                           jsize = jsize, jcolvec = jcolvec))
  multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
}

# why H3K27me3 look the same for different genes? One cell has huge weights on topic 23
# jgene <- "Hoxc13"

# do Sox cluster, all numbers
ref.mark <- "H3K4me1"
jgrep <- "Sox"
jsub <- subset(out.objs[[ref.mark]]$regions.annot, grepl(jgrep, SYMBOL)) %>% 
  rowwise() %>%
  mutate(suffix = as.numeric(strsplit(SYMBOL, jgrep)[[1]][[2]])) %>%  # get number from string
  arrange(suffix)
jgenes <- unique(jsub$SYMBOL)
for (jgene in jgenes){
  out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
  jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]])
  m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], 
                                                           show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                           jsize = jsize, jcolvec = jcolvec))
  multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
}

