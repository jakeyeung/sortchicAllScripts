# Jake Yeung
# Date of Creation: 2019-02-05
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/plot_marker_genes_all_4_marks.R
# Given list of genes, find best peak (relative to a mark) and plot bin across all 4 marks

jstart <- Sys.time()

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

jbin <- FALSE
if (jbin){
  Kstr <- "5_10_15_20_25"
} else {
  Kstr <- "5_15_25"
}
infs <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jmark, ".CountThres0.K-", Kstr, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
# out.objs <- lapply(jmarks, LoadLDABins)
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
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

# can save time by precalculating the UMAP and feeding it into the plot functions, also can customize for each UMAP



# regions.annot <- out.objs[[1]]$regions.annotated

outdir <- "~/Dropbox/scCHiC_figs/FIG4_BM/marker_genes"
dir.create(outdir)
fname <- "marker_genes.pdf"
pdf(file.path(outdir, fname), useDingbats = FALSE)


# find neutrophilmarkers
ref.mark <- "H3K4me1"
jgene <- "S100a8"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)

# Hbb gene
ref.mark <- "H3K4me1"
jgene <- "Hbb"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 50000)
jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
m.lst1 <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], use.count.mat = log(out.objs[[jmark]]$count.mat + 1),
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
m.lst2 <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], use.count.mat = NULL,
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst1[[1]], m.lst1[[2]], m.lst1[[3]], m.lst1[[4]], 
          m.lst2[[1]], m.lst2[[2]], m.lst2[[3]], m.lst2[[4]],
          cols = 4)

# Random gene
ref.mark <- "H3K4me1"
jgene <- "RandomlyPicked"
jpeak <- sample(out.objs[[ref.mark]]$regions.annot$region_coord, size = 1)
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], use.count.mat = out.objs[[jmark]]$count.mat, 
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], use.count.mat = NULL, 
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)


# Sox6 gene
# ref.mark <- "H3K4me1"
jgene <- "Sox6"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 100000)

# jpeak <- "chr7:115400000-115500000"
jpeak <- "chr7:115420000-115520000"
# jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)

# Car1
jgene <- "Car1"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 10000)
jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)


jgenes <- c("Cebpe", "Cebpa", "Elane", "Prtn3", "Mpo", "Flt3", "Ifitm1", "Lmo4", "Ccl5", "Prss34", "Meis1", "Cd74", "Gata2", "Car1", "Car2")

for (jgene in jgenes){
  out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 10000)
  jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
  m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                           show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                           jsize = jsize, jcolvec = jcolvec))
  multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
}

# Plot Hoxc13

ref.mark <- "H3K27me3"
jgene <- "Hoxc13"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
(jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]]))
m.lst1 <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], 
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
m.lst2 <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], use.count.mat = log(out.objs[[jmark]]$count.mat + 1),  
                                                         show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                         jsize = jsize, jcolvec = jcolvec))
multiplot(m.lst1[[1]], m.lst1[[2]], m.lst1[[3]], m.lst1[[4]], 
          m.lst2[[1]], m.lst2[[2]], m.lst2[[3]], m.lst2[[4]],
          cols = 4)

multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)

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
  print(jpeak)
  # m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]], 
                                                           # show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene,
                                                           # jsize = jsize, jcolvec = jcolvec))
  # multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
}

# why H3K27me3 look the same for different genes? One cell has huge weights on topic 23
# jgene <- "Hoxc13"



dev.off()

print(Sys.time() - jstart)