# Jake Yeung
# Date of Creation: 2020-10-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/H3K27me3_merged/cleanup_data_for_LDA.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(irlba)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jmarks <- c("H3K27me3" = "H3K27me3")
jmark <- jmarks[[1]]

inf.chromo <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs/countTablesAndRZr1only_ByChromo.NewFilters/PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230.sorted.tagged.countTable.ByChromo.WithSpikeIns.NoChromo.csv"
inf.bin <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs/countTablesAndRZr1only_ByChromo.NewFilters/PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230.sorted.tagged.countTable.binsize_50000.csv"
inf.rz <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs/countTablesAndRZr1only_ByChromo.NewFilters/PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230.sorted.tagged.countTable.RZ.csv"

min.total <- 1000
max.total <- 50000
min.l2r <- 0
max.l2r <- 5
min.TAfrac <- 0.6

jtitle.countsfilt <- paste("MinChromoCuts:", min.total, "MaxChromoCuts:", max.total)
jtitle.l2rfilt <- paste("MinLogRatio:", min.l2r, "MaxLogRatio:", max.l2r)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.merged_across_runs"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("qc_plots.", Sys.Date(), ".minl2r_", min.l2r, ".pdf"))
# assertthat::assert_that(file.exists(outpdf))

pdf(file = outpdf, useDingbats = FALSE)

# Get good cells ----------------------------------------------------------

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.rz <- ReadLH.SummarizeTA(inf.rz)
# dat.chromo <- ReadChrReads(inf.chromo)
dat.chromo <- GetChromoCounts(inf.chromo, spikeinchromo = jspikeinchromo, chromos.keep = jchromos) %>%
  filter(chromo == "1")
dat.bin <- ReadMatSlideWinFormat(inf.bin)

dat.rzchromo <- left_join(dat.rz, dat.chromo)

dat.rzchromo <- dat.rzchromo %>%
  rowwise() %>%
  mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
         colcoord = AddPlateCoordinates(samp)$colcoord,
         is.empty = rowcoord <= 8 & colcoord == 1,
         mark = jmark)

dat.rz.merged <- dat.rzchromo
# dat.rz.merged$mark <- jmark

dat.out.lst <- list()
dat.out.lst[[jmark]] <- list(dat.rz = dat.rzchromo, dat.mat = dat.bin)


ggplot(dat.rzchromo, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(ylim = c(0, 1))

ggplot(dat.rzchromo, aes(x = log10(chromocounts), y = TA.frac, color = is.empty)) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_cartesian(ylim = c(0, 1))



dat.rz.merged <- lapply(dat.out.lst, function(jdat){
  jdat$dat.rz
}) %>%
  bind_rows()

# dat.rz.merged <- dat.rz.merged %>%
#   rowwise() %>%
#   mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
#          colcoord = AddPlateCoordinates(samp)$colcoord)
# dat.rz.merged <- dat.rz.merged %>%
#   rowwise() %>%
#   mutate(is.empty = rowcoord <= 8 & colcoord == 1)


m1 <- ggplot(dat.rz.merged, aes(x = log10(chromocounts), y = TA.frac, color = is.empty, size = is.empty)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1)

m2 <- ggplot(dat.rz.merged, aes(x = log10(spikeincounts), y = TA.frac, color = is.empty, size = is.empty)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1)

m3 <- ggplot(dat.rz.merged, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, size = is.empty, color = is.empty)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1)

m4 <- ggplot(dat.rz.merged, aes(x = spikeincounts, y = chromocounts, size = is.empty, color = is.empty)) + 
  geom_point(alpha = 0.25)  + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1)

JFuncs::multiplot(m1, m2, m3, cols = 3)

print(m4)





# Filter two ways and run LDA  --------------------------------------------

# without spikeins


dat.rz.merged <- dat.rz.merged %>%
  rowwise() %>%
  mutate(is.good = chromocounts >= min.total & chromocounts <= max.total & TA.frac >= min.TAfrac & !is.empty,
         is.good.l2r = log2(chromocounts / spikeincounts) >= min.l2r & log2(chromocounts / spikeincounts) >= min.l2r & TA.frac >= min.TAfrac & !is.empty)


m1 <- ggplot(dat.rz.merged, aes(x = log10(chromocounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) +
  ggtitle(jtitle.countsfilt)

m2 <- ggplot(dat.rz.merged, aes(x = log10(spikeincounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.countsfilt)

m3 <- ggplot(dat.rz.merged, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.good)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.countsfilt)

m4 <- ggplot(dat.rz.merged, aes(x = spikeincounts, y = chromocounts, color = is.good)) + 
  geom_point(alpha = 0.25)  + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.countsfilt)

JFuncs::multiplot(m1, m2, m3, cols = 3)

print(m4)



m1 <- ggplot(dat.rz.merged, aes(x = log10(chromocounts), y = TA.frac, color = is.good.l2r)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.l2rfilt)

m2 <- ggplot(dat.rz.merged, aes(x = log10(spikeincounts), y = TA.frac, color = is.good.l2r)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.l2rfilt)

m3 <- ggplot(dat.rz.merged, aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.good.l2r)) + 
  geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.l2rfilt)

m4 <- ggplot(dat.rz.merged, aes(x = spikeincounts, y = chromocounts, color = is.good.l2r)) + 
  geom_point(alpha = 0.25)  + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(jtitle.l2rfilt)

JFuncs::multiplot(m1, m2, m3, cols = 3)

print(m4)



# Make cut tables ---------------------------------------------------------

cells.keep <- subset(dat.rz.merged, is.good)$samp
cells.keep.l2r <- subset(dat.rz.merged, is.good.l2r)$samp

mat.filt.lst <- lapply(dat.out.lst, function(jdat){
  jmat <- jdat$dat.mat
  cells.keep <- colnames(jmat) %in% cells.keep
  return(jmat[, cells.keep])
})

mat.filt.l2r.lst <- lapply(dat.out.lst, function(jdat){
  jmat <- jdat$dat.mat
  cells.keep <- colnames(jmat) %in% cells.keep.l2r
  return(jmat[, cells.keep])
})


# Check LSI ---------------------------------------------------------------


lapply(jmarks, function(jmark){
  print(jmark)
  lsi.out <- RunLSI(as.matrix(mat.filt.lst[[jmark]]))
  dat.lsi <- DoUmapAndLouvain(lsi.out$u, jsettings)
  m <- ggplot(dat.lsi, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) +
    ggtitle(paste("LSI is good:", jmark))
  print(m)
})


lapply(jmarks, function(jmark){
  print(jmark)
  lsi.out <- RunLSI(as.matrix(mat.filt.l2r.lst[[jmark]]))
  dat.lsi <- DoUmapAndLouvain(lsi.out$u, jsettings)
  m <- ggplot(dat.lsi, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) +
    ggtitle(paste("LSI l2r filt:", jmark))
  print(m)
})



# Write counts ------------------------------------------------------------

lapply(jmarks, function(jmark){
  print(jmark)
  jmat <- mat.filt.lst[[jmark]]
  outname <- paste0("count_mat_", jmark, "_counts_filt.", Sys.Date(), ".minl2r_", min.l2r, ".rds")
  outf <- file.path(outdir, outname)
  saveRDS(jmat, file = outf)
})

lapply(jmarks, function(jmark){
  print(jmark)
  jmat <- mat.filt.l2r.lst[[jmark]]
  print(dim(jmat))
  outname <- paste0("count_mat_", jmark, "_l2r_filt.", Sys.Date(), ".minl2r_", min.l2r, ".rds")
  outf <- file.path(outdir, outname)
  saveRDS(jmat, file = outf)
})


dev.off()