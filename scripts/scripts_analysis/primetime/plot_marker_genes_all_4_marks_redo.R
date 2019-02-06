# Jake Yeung
# Date of Creation: 2019-02-05
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/plot_marker_genes_all_4_marks_redo.R
# Redo with norm counts

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


# Functions ---------------------------------------------------------------

PlotUmapAllMarks <- function(jmarks, tm.result.lst, jpeak, juse.count.mat, dat.umap.lst, jgene, jsize, jcolvec, .log=TRUE){
  m.lst <- lapply(jmarks, function(jmark) PlotImputedPeaks2(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                            use.count.mat = NULL,
                                                            usettings=dat.umap.lst[[jmark]], 
                                                            gname = jgene,
                                                            jsize = jsize, jcolvec = jcolvec, .log = TRUE))
  m.lst2 <- lapply(jmarks, function(jmark) PlotImputedPeaks2(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],
                                                             use.count.mat = juse.count.mat[[jmark]],
                                                             usettings=dat.umap.lst[[jmark]], 
                                                             gname = jgene,
                                                             jsize = jsize, jcolvec = jcolvec, .log = TRUE))
  multiplot(m.lst[[1]], m.lst2[[1]], 
            m.lst[[2]], m.lst2[[2]], 
            m.lst[[3]], m.lst2[[3]],
            m.lst[[4]], m.lst2[[4]],
            cols = 4)
}

# load from bins ----------------------------------------------------------

jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
jcolvec <- c("blue", "white", "red")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

jbin <- TRUE
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

count.mat.lst <- lapply(out.objs, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

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
# custom settings for each UMAP
nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.4, 0.2, 0.2)
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed), 
                              nn.vec, jmindist.vec, SIMPLIFY = FALSE)
# plot umaps to check
topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks

umap.plots <- lapply(dat.umap.lst, function(dat.umap){
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2])
  m <- ggplot(dat, aes(x = umap1, y = umap2)) + 
    geom_point(size = 0.5, alpha = 0.2) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})
multiplot(umap.plots[[1]], umap.plots[[2]], umap.plots[[3]], umap.plots[[4]], cols = 4)

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
PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)

# Hbb gene
ref.mark <- "H3K4me1"
jgene <- "Hbb"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 50000)
jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)

# Random gene
set.seed(jseed)
ref.mark <- "H3K4me1"
gene <- "RandomlyPickedPeak"
jpeak <- sample(out.objs[[ref.mark]]$regions.annot$region_coord, size = 1)
PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)


# Sox6 gene
# ref.mark <- "H3K4me1"
# jgene <- "Sox6"
# out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 100000)
# 
# jpeak <- "chr7:115400000-115500000"
jgene <- "Sox6"
jpeak <- "chr7:115420000-115520000"  # hand-picked by Alexander
# jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)

# Car1
jgene <- "Car1"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 10000)
jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)


jgenes <- c("Cebpe", "Cebpa", "Elane", "Prtn3", "Mpo", "Flt3", "Ifitm1", "Lmo4", "Ccl5", "Prss34", "Meis1", "Cd74", "Gata2", "Car1", "Car2")

for (jgene in jgenes){
  out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot, dist = 10000)
  jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])
  PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)
}

# Plot Hoxc13

ref.mark <- "H3K27me3"
jgene <- "Hoxc13"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
(jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]]))
PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)

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
  PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)
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
  PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE)
}
dev.off()

print(Sys.time() - jstart)