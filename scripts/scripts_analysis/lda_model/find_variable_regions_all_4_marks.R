# Jake Yeung
# Date of Creation: 2019-01-21
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/find_variable_regions_all_4_marks.R
# Get variable regions for all 4 marks

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


source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

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


jchips <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
out.objs <- lapply(jchips, LoadLDA)
names(out.objs) <- jchips



# Marks -------------------------------------------------------------------

GetVar <- function(tm.result, regions.annotated){
  mat.norm <- t(tm.result$topics %*% tm.result$terms)
  # find peaks with largest range
  maxmin <- sort(apply(mat.norm, 1, function(jcol) mad(jcol)), decreasing = TRUE)
  # what are the sizes?
  peak.size <- sapply(names(maxmin), function(x) ParseCoord(x)$end - ParseCoord(x)$start)
  # normalize maxmin by peaksize
  maxmin.norm <- maxmin / peak.size
  dat.var <- data.frame(Var = maxmin, peak.size = peak.size, coord = names(maxmin))
  regions.annotated <- dplyr::left_join(regions.annotated, dat.var, by = c("region_coord"="coord"))
  return(regions.annotated)
}

regions.with.var <- lapply(out.objs, function(out){
  out$regions.annotated <- out$regions.annotated %>%
    mutate(region_coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"))
  return(GetVar(out$tm.result, out$regions.annotated))
})

# label assay and merge
for (m in names(regions.with.var)){
  regions.with.var[[m]]$mark <- m
}

regions.with.var.long <- bind_rows(regions.with.var)


# Find covarying regions? -------------------------------------------------


jsub <- subset(regions.with.var.long, abs(distanceToTSS) < 10000 & grepl(jgene, SYMBOL)) %>%
  group_by(mark) %>%
  filter(Var == max(Var)) %>%
  group_by(SYMBOL) %>%
  summarise(Var.mean = mean(Var),
            Var.mad = mad(Var))

jsub.var <- subset(regions.with.var.long, abs(distanceToTSS) < 10000) %>%
  group_by(SYMBOL, mark) %>%
  filter(Var == max(Var)) %>%
  group_by(SYMBOL) %>%
  filter(length(Var) == 4) %>%
  summarise(Var.mean = median(Var),
            Var.mad = mad(Var))

# find regions that have high mean and low variability/entropy


jgene <- "Actr3"
jgene <- "Herc4"
jgene <- "Nfia"
jgene <- "Galnt1"
jgene <- "Glg1"
jgenes <- c("Actr3", "Herc4", "Nfia", "Galnt1", "Glg1")
jdist <- 50000
# /pdf("~/Dropbox/scCHiC_figs/FIG4_BM/LDA_outputs/variable_regions_genes.pdf")
for (jgene in jgenes){
  system.time(
    plots.out <- PlotAllMarks(jgene, jchips, jdist, out.objs)
  )
  multiplot(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], cols = 4)
}
#
