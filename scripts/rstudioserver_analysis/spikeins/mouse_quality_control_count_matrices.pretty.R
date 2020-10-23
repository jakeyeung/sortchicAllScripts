# Jake Yeung
# Date of Creation: 2020-08-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse_quality_control_count_matrices.pretty.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(irlba)

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda"
# assertthat::assert_that(dir.exists(outdir))
# outf <- file.path(outdir, "H3K4me3_BM.rds")
# outf.filt <- file.path(outdir, "H3K4me3_BM.dev_filt.rds")
# outspikeins <- file.path(outdir, "H3K4me3_BM.spikeins.txt")
# outpdf <- file.path(outdir, "H3K4me3_BM.pdf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt"
dir.create(outdir)
outf1 <- file.path(outdir, "H3K4me3_BM.rds")
outpdf <- file.path(outdir, "H3K4me3_BM.pdf")
outf.devfilt <- file.path(outdir, "H3K4me3_BM.dev_filt.rds")
outf.spikeins <- file.path(outdir, "H3K4me3_BM.spikeins.txt")
outf.forproj <- file.path(outdir, "H3K4me3_BM.match_rownames_with_old.rds")

assertthat::assert_that(!file.exists(outf1))
assertthat::assert_that(!file.exists(outpdf))
assertthat::assert_that(!file.exists(outf.devfilt))
assertthat::assert_that(!file.exists(outf.spikeins))

pdf(outpdf, useDingbats = FALSE)

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
log2filt <- 0.1
TAfracmin <- 0.6



# Load spikein info -------------------------------------------------------

inf.bychromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.ByChromo.WithSpikeIns.csv")
jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")
dat.chromo <- GetChromoCounts(inf.bychromo, spikeinchromo = jspikeinchromo, chromos.keep = NA)

jdat <- subset(dat.chromo, chromo == 1) %>%
  dplyr::select(-chromo)
jdat <- as.data.frame(jdat)
rownames(jdat) <- jdat$samp



# Merge with TAfrac -------------------------------------------------------

dat.rz.merged <- left_join(dat.rz, jdat) %>%
  rowwise() %>%
  mutate(s2n = log2(chromocounts / spikeincounts),
         is.good = total.count >= cutsmin & TA.frac >= TAfracmin & s2n >= log2filt)

dat.rz.filt <- as.data.frame(subset(dat.rz.merged, is.good))
rownames(dat.rz.filt) <- dat.rz.filt$samp

cells.good <- dat.rz.filt$samp

# ggplot(dat.rz.merged, aes(x = log10(total.count), y = TA.frac, color = is.good)) + 
#   geom_point() + 
#   coord_cartesian(ylim = c(0, 1)) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.merged, aes(x = log10(chromocounts), y = TA.frac, color = is.good)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0, 1)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.merged, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.good)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(ylim = c(0, 1)) + 
  geom_vline(xintercept = 0.1, linetype = "dotted", color = "grey55")



# Load mats ---------------------------------------------------------------

inf.mat <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.binsize_50000.csv")

mat <- ReadMatSlideWinFormat(inf.mat)

cols.i <- colnames(mat) %in% cells.good

print(dim(mat))
mat.filt <- mat[, cols.i]
print(dim(mat.filt))

lsi.out <- RunLSI(as.matrix(mat.filt))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(lsi.out$u, jsettings)

m2 <- ggplot(dat.umap, aes(x = umap1, y = umap2)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m2)


# Deviance filter ---------------------------------------------------------

jdat.reordered <- dat.rz.filt[colnames(mat.filt), ]
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

dim(mat.filt)
dim(mat.filt.gdev)

dev.off()


# Filter for projections --------------------------------------------------

inf.ref <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.for_projections/H3K4me3_padded_zeros_for_projections.rds"

mat.ref <- readRDS(inf.ref)

rows.filt1 <- rownames(mat.filt) %in% rownames(mat.ref)

mat.filt.forproj1 <- mat.filt[rows.filt1, ]

# pad zeros 
rnames.add <- rownames(mat.ref)[!rownames(mat.ref) %in% rownames(mat.filt)]

mat.add <- matrix(data = 0, nrow = length(rnames.add), ncol = ncol(mat.filt), dimnames = list(rnames.add, colnames(mat.filt)))

mat.filt.forproj <- Matrix(rbind(mat.filt.forproj1, mat.add), sparse = TRUE)

# reearrange rows
mat.filt.forproj <- mat.filt.forproj[rownames(mat.ref), ]


# Write new count tables  -------------------------------------------------


saveRDS(mat.filt, file = outf1)
saveRDS(mat.filt.gdev, file = outf.devfilt)
saveRDS(mat.filt.forproj, file = outf.forproj)
fwrite(dat.rz.filt, file = outf.spikeins)






