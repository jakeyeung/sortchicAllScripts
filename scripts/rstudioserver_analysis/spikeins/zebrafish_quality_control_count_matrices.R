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

species <- "zebrafish"
outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_", species, "_for_lda")
dir.create(outdir)
assertthat::assert_that(dir.exists(outdir))
outf <- file.path(outdir, "H3K4me3_WKM.rds")
outpdf <- file.path(outdir, "H3K4me3_WKM.pdf")

# Load mats ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


inf.rz <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/VAN4969/", species, "/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/ZF-WKM-EtOH-H3K4me3-2.sorted.tagged.countTable.RZ.csv"))

# inf.rz <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/VAN4969/", species, "/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.RZ.csv"))
inf.mat <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/VAN4969/", species, "/tagged_bams/countTablesAndRZr1only.NewFilters/ZF-WKM-EtOH-H3K4me3-2.sorted.tagged.countTable.binsize_50000.csv"))

dat.rz <- ReadLH.SummarizeTA(inf.rz)

cutsmin <- 500
TAfracmin <- 0.6

m1 <- ggplot(dat.rz, aes(x = total.count, y = TA.frac)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + 
  geom_vline(xintercept = cutsmin) + 
  geom_hline(yintercept = TAfracmin)

print(m1)


dat.rz.filt <- subset(dat.rz, total.count >= cutsmin & TA.frac >= TAfracmin)

cells.good <- dat.rz.filt$samp

# Load mats ---------------------------------------------------------------


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

pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
  print(m1)
  print(m2)
dev.off()

# Write table  ------------------------------------------------------------

saveRDS(mat.filt, file = outf)



