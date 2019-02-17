# Jake Yeung
# Date of Creation: 2019-02-08
# File: ~/projects/scchic/scripts/scripts_analysis/mara_downstream/mara_downstream.R
# Look at MARA output

library(scales)
library(data.table)
library(reshape2)

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
library(forcats)
library(ggrepel)
library(biomaRt)

library(igraph)  # louvain

library(Gviz)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

source("scripts/Rfunctions/PlotFunctions.R")

# Constants you can tweak -------------------------------------------------



jchip <- "H3K4me1"





# settings for UMAP
if (jchip == "H3K4me1"){
  nn=40
} else {
  nn=35
}
# nn=15
nnterms <- 15
jmetric='euclidean' 
jmindist=0.4
jseed=123

plotout <- paste0("/tmp/", jchip, "_LDA_bins_top_regions.pdf")

# LDA was run on binarized matrix or not. 
# I was thinking this binarized matrix would help reduce weird genomic regions with way too many reads. 
# Because we expect the count matrix to have only a few reads per bin per cell. 
# Can tweak this to TRUE or FALSE
jbin <- "TRUE"



# Load MARA ---------------------------------------------------------------



# mdir <- "/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me1.filt_0.99/hiddenDomains_motevo_merged.closest.long.scale_1.center_0/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me1.filt_0.99"
mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
               jchip, 
               ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
               "/",
               "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jchip, ".filt_0.99.center_TRUE")
assertthat::assert_that(dir.exists(mdir))



act.mat <- fread(file.path(mdir, "Activities"), header = FALSE)
se.mat <- fread(file.path(mdir, "StandardError"), header = FALSE)
cnames <- unlist(fread(file.path(mdir, "Colnames"), header = FALSE), use.names = FALSE)
zscores <- fread(file.path(mdir, "Zscores"), header = FALSE)

colnames(zscores) <- c("motif", "zscore")
zscores <- zscores %>% arrange(desc(zscore))


# Load --------------------------------------------------------------------

dirmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"


if (jbin){
  inf <- file.path(dirmain, paste0("lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_10_15_20_25.Robj"))
} else {
  inf <- file.path(dirmain, paste0("lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_15_25.Robj"))
}

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

out.lda <- ChooseBestLDA(out.lda)
(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)

topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# check your umap settings
jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jchip, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)

print(sessionInfo())

# Plot dat umap -----------------------------------------------------------
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75)
}, jcol.rgbs, seq(kchoose))





# Combine the two ---------------------------------------------------------





# plot motif activity 
dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))

barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))

cnames.new <- sapply(cnames, function(x){
  jmark <- strsplit(x, "\\.")[[1]][[1]]
  jtiss <- strsplit(x, "\\.")[[1]][[2]]
  jmouse <- strsplit(x, "\\.")[[1]][[3]]
  jrep <- paste0("rep", strsplit(x, "\\.")[[1]][[4]])
  jcell <- cellhash.bc[[strsplit(x, "\\.")[[1]][[6]]]]
  xnew <- paste(jtiss, jmark, jmouse, jrep, jcell, sep = "_")
})

colnames(act.mat) <- c("motif", cnames.new)
colnames(se.mat) <- c("motif", cnames.new)

act.long <- tidyr::gather(act.mat, -motif, key = "cell", value = "activity")

# dat.merged <- left_join(act.long, dat.umap.long)  # act.long has some bad cells??
dat.merged <- left_join(dat.umap.long, act.long)


jmotifs <- zscores$motif[1:40]
jmotifs <- c(jmotifs, "Tal1")

plotf <- paste0("~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/mara/", jchip, "_motifs_init.pdf")
pdf(plotf, useDingbats = FALSE)

for (jmotif in jmotifs){
  # jmotif <- zscores$motif[[15]]
  dat.sub <- subset(dat.merged, motif == jmotif)
  jzscore <- signif(subset(zscores, motif == jmotif)$zscore, digits = 2)
  jtitle <- paste(jmotif, "zscore:", jzscore)
  # dat.sub <- subset(dat.merged, motif == "Tal1")
  m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = activity)) + geom_point() + 
    ggtitle(jtitle) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = mean(dat.sub$activity))
  print(m)
}

dev.off()