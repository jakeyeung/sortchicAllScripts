# Jake Yeung
# Date of Creation: 2019-04-05
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/plot_umap_with_total_counts.R
# Plot umap with total counts

rm(list=ls())

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

source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")

# Functions ---------------------------------------------------------------

PlotUmapWithCellSums <- function(jmark, dat.umap.long.new.lst, count.mat.lst, jscale = 1){
  cell.sums <- data.frame(cell = colnames(count.mat.lst[[jmark]]), cellsum = Matrix::colSums(count.mat.lst[[jmark]]))
  # add it to umap and plot
  dat.umap <- left_join(dat.umap.long.new.lst[[jmark]], cell.sums)
  head(dat.umap.long.new.lst[[jmark]])
  cell.sums <- data.frame(cell = colnames(count.mat.lst[[jmark]]), cellsum = Matrix::colSums(count.mat.lst[[jmark]]))
  # add it to umap and plot
  dat.umap <- left_join(dat.umap.long.new.lst[[jmark]], cell.sums) 
  m1 <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = cellsum)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m2 <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = log10(cellsum))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  print(m2)
  return(dat.umap)
}

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"

load(inf, v=T)



# Get counts total and append to UMAP  ------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
jscales <- c(-1, -1, 1, 1)

jmark <- jmarks[[1]]
jmark <- jmarks[[2]]
jmark <- jmarks[[3]]
jmark <- jmarks[[4]]
jscale <- 1

# Load count.mat again unnormalized ---------------------------------------

infs.nobin <- list("/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K4me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-20_25_30_35.Robj",
                   "/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K9me3.CountThres0.K-20_25_30_35.Robj")
names(infs.nobin) <- c(jmarks)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
count.mat.lst <- lapply(out.objs.nobin, function(x) x$count.mat)

head(dat.umap.long.new.lst[[jmark]])
cell.sums <- data.frame(cell = colnames(count.mat.lst[[jmark]]), cellsum = Matrix::colSums(count.mat.lst[[jmark]]))
# add it to umap and plot
dat.umap <- left_join(dat.umap.long.new.lst[[jmark]], cell.sums)

pdf(paste0("~/data/scchic/pdfs/umap_with_counts.", Sys.Date(), ".pdf"), useDingbats = FALSE)
outs <- mapply(function(jmark, jscale) PlotUmapWithCellSums(jmark, dat.umap.long.new.lst, count.mat.lst, jscale), jmarks, jscales, 
               MoreArgs = list(dat.umap.long.new.lst = dat.umap.long.new.lst, count.mat.lst = count.mat.lst))
dev.off()
