# Jake Yeung
# Date of Creation: 2019-03-14
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/plot_marker_genes_repressive_marks_lessmem.R
# Redo with less memory??


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
library(scales)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")


jscale.fac <- 10^6
jpseudo <- 1

# Functions ---------------------------------------------------------------



# load from bins ----------------------------------------------------------

jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
# jcolvec <- c("blue", "gray80", "red")
jcolvec <- c("gray95", "gray50", "darkblue")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

jmarks.keep <- c("H3K27me3", "H3K9me3")
names(jmarks.keep) <- jmarks.keep

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

# infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
infs <- infs.nobin[jmarks.keep]
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
names(out.objs) <- jmarks
names(out.objs.nobin) <- jmarks

tm.result.lst <- lapply(out.objs[jmarks.keep], function(x) posterior(x[['out.lda']]))

# use nobin for mat
count.mat.lst <- lapply(out.objs.nobin[jmarks.keep], function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

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
nn.vec <- c(40, 40)
jmindist.vec <- c(0.2, 0.2)
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed), 
                              nn.vec, jmindist.vec, SIMPLIFY = FALSE)
# plot umaps to check
topics.mat.lst <- lapply(out.objs[jmarks.keep], function(x) x$tm.result$topics)
dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks.keep

umap.plots <- lapply(dat.umap.lst, function(dat.umap){
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2])
  m <- ggplot(dat, aes(x = umap1, y = umap2)) + 
    geom_point(size = 0.5, alpha = 0.2) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

# top hits 
top.peaks <- tidytext::tidy(out.objs[[2]]$out.lda, matrix = "beta") %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

# Plot all 4 marks for the shape 




# regions.annot <- out.objs[[1]]$regions.annotated

outdir <- "~/Dropbox/scCHiC_figs/FIG4_BM/marker_genes"
dir.create(outdir)
fname <- "marker_genes_redo_nofilt_rescale_repressive_lessmem.pdf"
# pdf(file.path(outdir, fname), useDingbats = FALSE)

multiplot(umap.plots[[1]], umap.plots[[2]], cols = length(infs))

# find neutrophilmarkers

# Random gene
set.seed(jseed)
ref.mark <- "H3K9me3"
jgene <- "IgH region"
jpeak <- "chr12:115560000-115660000"
PlotUmapAllMarks(jmarks.keep, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)


ref.mark <- "H3K27me3"
jgene <- "Hoxc13"
out.sub <- GetPeaksFromGene(jgene, out.objs[[ref.mark]]$regions.annot)
(jpeak <- SelectBestPeak(out.sub$peaks, out.objs[[ref.mark]]$regions.annot, tm.result.lst[[ref.mark]]))
PlotUmapAllMarks(jmarks.keep, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)

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
  PlotUmapAllMarks(jmarks, tm.result.lst, jpeak, juse.count.mat = count.mat.lst, dat.umap.lst, jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)
}

# dev.off()

print(Sys.time() - jstart)
