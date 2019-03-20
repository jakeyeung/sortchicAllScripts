# Jake Yeung
# Date of Creation: 2019-03-18
# File: ~/projects/scchic/scripts/scripts_analysis/make_primetime_objs/make_primetime_objs.R
# Make primetime objects, sick of loading every time 


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

library(igraph)

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
jcolvec <- c("gray90", "gray50", "darkblue")


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

# use nobin for mat
count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

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
nn.vec <- c(40, 35, 27, 40)
jmindist.vec <- c(0.2, 0.1, 0.2, 0.1)
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed), 
                              nn.vec, jmindist.vec, SIMPLIFY = FALSE)
# plot umaps to check
topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks

dat.umap.long.lst <- lapply(dat.umap.lst, function(dat.umap){
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
})


# Do louvains -------------------------------------------------------------

nn.louv.vec <- c(27, 27, 60, 100)
names(nn.louv.vec) <- jmarks
jmetric.louv='euclidean' 
jmindist.louv=0.4
jseed.louv=123
custom.settings.louv.lst <- lapply(nn.louv.vec, function(nn.louv) {
  jtmp <- GetUmapSettings(nn=nn.louv, 
                          jmetric=jmetric.louv, 
                          jmindist=jmindist.louv, 
                          seed = jseed.louv)
  return(jtmp)
})

dat.umap.long.lst <- mapply(DoLouvain, topics.mat.lst, custom.settings.louv.lst, dat.umap.long.lst, SIMPLIFY = FALSE)

# plot louvains with colors
jmark <- "H3K9me3"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
m.louvain <- ggplot(dat.umap.long.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
print(m.louvain)


# Write things to output --------------------------------------------------

outdir <- "~/data/scchic/robjs/primetime_objs"
dir.create(outdir)
save(dat.umap.long.lst, tm.result.lst, count.mat.lst, custom.settings.lst, custom.settings.louv.lst, topics.mat.lst, out.objs, out.objs.nobin, file = file.path(outdir, "four_marks_lda_output.RData"))

