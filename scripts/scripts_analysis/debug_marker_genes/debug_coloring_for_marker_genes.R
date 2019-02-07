# Jake Yeung
# Date of Creation: 2019-02-05
# File: ~/projects/scchic/scripts/scripts_analysis/debug_marker_genes/debug_coloring_for_marker_genes.R
# Consolidate why a random peak looks good and H3K27me3 sometimes show inverse between imputed and raw data.


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

# settings for UMAP
# single setting for all 4 
nn=40
nnterms <- 15
jmetric='euclidean'
jmindist=0.2
jseed=123
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)

# custom settings for each UMAP
nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.4, 0.2, 0.2)
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed), 
                              nn.vec, jmindist.vec)
# plot umaps to check
dat.umap <- umap(topics.mat, config = custom.settings)


# Plot Hox clusters for H3K27me3 ------------------------------------------

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

# multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 4)
multiplot(m.lst1[[1]], m.lst1[[2]], m.lst1[[3]], m.lst1[[4]], 
          m.lst2[[1]], m.lst2[[2]], m.lst2[[3]], m.lst2[[4]],
          cols = 4)

# check what the mat.norm range is for each plotimputedpeaks

# jmark <- "H3K27me3"
jmark.ref <- "H3K4me1"

tm.result <- tm.result.lst[[jmark.ref]]
jcounts.norm.lst <- lapply(jmarks, function(jmark){
  tm.result <- tm.result.lst[[jmark]]
  mat.norm <- t(tm.result$topics %*% tm.result$terms)
  row.i <- which(rownames(mat.norm) %in% jpeak)
  jcounts.norm <- mat.norm[row.i, ]
  # zscore?
  jcounts.norm <- scale(jcounts.norm, center = TRUE, scale = TRUE)
  return(jcounts.norm)
})

par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
for (jmark in jmarks){
  plot(density(jcounts.norm.lst[[jmark]]), main = jmark)
}
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

# plot imputation with zscore?

# change color scheme here
jcounts.norm <- jcounts.norm.lst[[jmark.ref]]
# jcol.counts <- ColorsByCounts(scale(jcounts.norm.lst[[3]]), nbreaks = 100, colvec = jcolvec)
# plot topics soft clustering weights
topics.mat <- tm.result$topics

dat.umap <- umap(topics.mat, config = custom.settings)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# jpeaks.str <- paste(peaks.keep.kb, collapse = ",")
# jmain <- paste0(jchip, " ", jlab, ' npeaks ', length(row.i), "\n", gname)
# prepare plot object
dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                  umap2 = dat.umap$layout[, 2], 
                  counts.norm = jcounts.norm)
m <- ggplot(dat, aes(x = umap1, y = umap2, col = jcounts.norm)) + 
  geom_point(size = jsize) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size=7), legend.position = "bottom") + 
  scale_color_gradientn(colours = c("blue", "white", "red"), values = scales::rescale(c(-2, 1, 8)), space = "Lab")
print(m)
  scale_color_gradient2(low = "blue", mid = "white", high = "red", space="Lab", limits = c(-2, 8))


