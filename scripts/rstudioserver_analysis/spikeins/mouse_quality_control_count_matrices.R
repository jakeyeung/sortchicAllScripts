# Jake Yeung
# Date of Creation: 2020-08-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse_quality_control_count_matrices.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(irlba)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda"
assertthat::assert_that(dir.exists(outdir))
outf <- file.path(outdir, "H3K4me3_BM.rds")
outf.filt <- file.path(outdir, "H3K4me3_BM.dev_filt.rds")
outspikeins <- file.path(outdir, "H3K4me3_BM.spikeins.txt")
outpdf <- file.path(outdir, "H3K4me3_BM.pdf")

# Load mats ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.rz <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.RZ.csv")

dat.rz <- ReadLH.SummarizeTA(inf.rz)

m1 <- ggplot(dat.rz, aes(x = total.count, y = TA.frac)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10()

print(m1)

cutsmin <- 1000
TAfracmin <- 0.6

dat.rz.filt <- subset(dat.rz, total.count >= cutsmin & TA.frac >= TAfracmin)

cells.good <- dat.rz.filt$samp

# Load mats ---------------------------------------------------------------

inf.mat <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.binsize_50000.csv")

mat <- ReadMatSlideWinFormat(inf.mat)

cols.i <- colnames(mat) %in% cells.good
mat.filt <- mat[, cols.i]

lsi.out <- RunLSI(as.matrix(mat.filt))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(lsi.out$u, jsettings)

m2 <- ggplot(dat.umap, aes(x = umap1, y = umap2)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

if (!file.exists(outpdf)){
  pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
    print(m1)
    print(m2)
  dev.off()
}

# Write table  ------------------------------------------------------------

if (!file.exists(outf)){
  saveRDS(mat.filt, file = outf)
}


# Write spikeincounts -----------------------------------------------------

# inf.bychromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.ByChromo.csv")
inf.bychromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.ByChromo.WithSpikeIns.csv")

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.filt.long <- GetChromoCounts(inf.bychromo, spikeinchromo = jspikeinchromo, chromos.keep = NA)

jdat <- subset(dat.filt.long, chromo == 1 & samp %in% cells.good) %>%
  dplyr::select(-chromo)
jdat <- as.data.frame(jdat)
rownames(jdat) <- jdat$samp

# dat.spikeins.mat <- data.frame(samp = dat.filt.long )

# summarize chromo to spkeincounts per cell 
plot(density(log10(jdat$chromocounts)))
plot(density(log10(jdat$spikeincounts)))
plot(density(log2(jdat$chromocounts / jdat$spikeincounts)))



# save to output as txt
if (!file.exists(outspikeins)){
  fwrite(jdat, file = outspikeins)
} else {
  print(paste(outspikeins, "exists... skipping writing outspikeins"))
}


# top5000 -----------------------------------------------------------------

jdat.reordered <- jdat[colnames(mat.filt), ]

nvec <- jdat.reordered$spikeincounts

gdevs <- apply(mat.filt, 1, function(xvec){
  scchicFuncs::binomial_deviance(x = xvec, p = sum(xvec) / sum(nvec), n = nvec)
})

plot(density(gdevs))
plot(density(log10(gdevs)))


# check hihg dev gene
topg.i <- which.max(gdevs)
topg <- rownames(mat)[topg.i]
cdat <- data.frame(count = mat[topg, ], cell = colnames(mat), stringsAsFactors = FALSE)

plot(density(cdat$count), main = topg, xlab = "counts")

topn <- 5000
genes.keep.vec <- sort(gdevs, decreasing = TRUE)[1:topn]
genes.keep <- names(genes.keep.vec)

range(gdevs[genes.keep])

mat.filt.gdev <- mat.filt[genes.keep, ]

if (!file.exists(outf.filt)){
  saveRDS(mat.filt.gdev, file = outf.filt)
}

