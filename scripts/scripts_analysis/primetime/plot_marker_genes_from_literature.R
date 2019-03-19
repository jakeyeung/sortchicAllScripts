# Jake Yeung
# Date of Creation: 2019-03-18
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/plot_marker_genes_from_literature.R
# Find some literature genes, see what they look like in our dataset 


rm(list=ls())

setwd("~/projects/scchic")

jstart <- Sys.time()


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

library(cowplot)

library(ggrastr)

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

library(Gviz)

library(biomaRt)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")

source("scripts/Rfunctions/IntegrateData.R")


# Load active marks LDA output --------------------------------------------


# marks.keep <- c("H3K4me1", "H3K4me3")
# marks.keep <- c("H3K27me3", "H3K9me3")
marks.keep <- c("H3K4me1", "H3K27me3")
jmark <- marks.keep[[1]]  # reference mark for later 

jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
# jcolvec <- c("blue", "gray80", "red")
jcolvec <- c("gray70", "gray50", "darkblue")


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

mat.impute.lst <- lapply(tm.result.lst, function(tm.result) t(tm.result$topic %*% tm.result$term))

# use nobin for mat
count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

print(lapply(mat.impute.lst, dim))

# Wrangle data ------------------------------------------------------------

# top.thres <- 0.999
top.regions.lst <- lapply(jmarks, function(jmark){
  topic.regions <- out.objs[[jmark]]$topic.regions  # for each cluster
  # topic.regions <- lapply(out.objs[[jmark]]$out.lda, function(clst){
  #   return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  # })
  top.regions <- unique(unlist(topic.regions))
  
  top.regions <- unique(unlist(topic.regions))
  return(top.regions) 
})
top.regions.merge <- unique(unlist(top.regions.lst))

mat.merged.lst <- lapply(mat.impute.lst, function(mat.impute){
  row.i <- which(top.regions.merge %in% rownames(mat.impute))
  return(as.data.frame(mat.impute[row.i, ]))
})

rnames.lst <- lapply(mat.merged.lst, function(x) rownames(x))
top.regions.intersect <- purrr::reduce(.x = rnames.lst, .f = intersect)

print(length(top.regions.intersect))

mat.merged.lst <- lapply(mat.merged.lst, function(mat.merged){
  return(mat.merged[top.regions.intersect, ])
})


# Find some genes ---------------------------------------------------------

jscale.fac <- 1
jpseudo <- 10^-6

# can save time by precalculating the UMAP and feeding it into the plot functions, also can customize for each UMAP
# custom settings for each UMAP
jmetric='euclidean'
jmindist=0.2
jseed=123
nn.vec <- c(40, 35, 40, 40)
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

jpeak <- "chr7:73520000-73620000"  # top loads1 peak

# Plot some interesting genes  --------------------------------------------

# correlate cells to cells

# find S100a8 genes and plot repressive marks 

jpeak <- "chr3:90660000-90760000"

jmark <- "H3K4me1"
marks.keep <- c("H3K4me1", "H3K27me3")
print(jpeak)
# jpeak <- strsplit(jpeak, split = ";")[[1]][[1]]
# (jgene <- subset(out.objs[[jmark]]$regions.annot, region_coord == jpeak)$SYMBOL[[1]])
jgene <- "Tal1"

jgene <- "Rprl1"
jgene <- "S100a7a"

jgene <- "Hoxc13"

jgene <- "Ezh2"

jgene <- "Ly6c1"
jgene <- "Ly6c2"
jgene <- "Ankrd28"

jgene <- "Ulk2"

jgene <- "Ncor1"

annot <- GetPeaksFromGene(jgene, out.objs[[jmark]]$regions.annot, dist = 10^5)
(jpeak <- annot$peaks[[1]])

system.time(
  PlotUmapAllMarks(jmarks[marks.keep], tm.result.lst[marks.keep], jpeak, juse.count.mat = count.mat.lst[marks.keep], dat.umap.lst[marks.keep], jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)
)
# 
# # Plot GTracks for region of interest
# 
# x.sum <- readRDS("~/data/scchic/robjs/x_sum.H3K9me3.rds")
# 
# mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
# jstart <- as.numeric(GetStart(jpeak)) - 6*10^5
# jend <- as.numeric(GetEnd(jpeak)) + 6*10^5
# jchromo <- GetChromo(jpeak)
# 
# system.time(
#   PlotGTrack2(x.sum, jstart, jend, mart.obj, louvains, gen = "mm10", chr = jchromo, jheight = "auto")
# ) 
# 
# 
# 

