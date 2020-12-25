# Jake Yeung
# Date of Creation: 2020-10-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/filter_good_cells_for_bam_clustering.G1only.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)


# Load mats ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

bsize <- "50000"
jthres <- 1.5
# jthres <- 2
jsuffix <- paste0("nonzerofracthres_", jthres, ".bsize_", bsize)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K562_clean_by_G1filt.fracnonzerofilt"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("badcells_filt_by_nonzerofrac.", Sys.Date(), ".", jsuffix, ".pdf"))

pdf(outpdf, useDingbats = FALSE)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2.G1filt"

indir.chromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_TAfrac.NewFilters")


mats <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("K562_count_tables_", bsize, ".", jmark, ".G1filt.rds")
  inf <- file.path(indir, fname)
  readRDS(inf)
})


dats.meta <- lapply(jmarks, function(jmark){
  mat <- mats[[jmark]]
  dat.meta <- data.frame(samp = colnames(mat), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(plate = ClipLast(samp, jsep = "_"),
           mark = jmark)
  return(dat.meta)
})
  
# Add meta data -----------------------------------------------------------

dats.meta.annot <- lapply(dats.meta, function(dat){
  dat <- AddCellCycleLabel.bydat(dat)
  dat$cellcycle.str <- "G1filt"
  return(dat)
})


# Write output ------------------------------------------------------------

dats.meta.annot.out <- lapply(dats.meta.annot, function(dat){
  dat <- dat %>%
    mutate(cluster = cellcycle.str,
           cell = samp) %>%
    dplyr::select(c(cell, cluster, plate, mark, rowcoord, colcoord))
})


# Add spikeins ------------------------------------------------------------




# Load chromos ------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

infs.chromo <- list.files(indir.chromo, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
assertthat::assert_that(length(infs.chromo) > 0)


jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos <- lapply(infs.chromo, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = NA)
}) %>%
  bind_rows()


# Load LH counts ----------------------------------------------------------

infs.lh <- list.files(indir.lh, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
dat.lh <- lapply(infs.lh, function(inf){
  dat.filt.long.lh <- ReadLH.SummarizeTA(inf)
}) %>%
  bind_rows() %>%
  mutate(experi = ClipLast(samp, jsep = "_"))

chromocounts <- subset(dat.chromos, chromo == "1", select = c(samp, chromocounts, spikeincounts))

dat.lh <- left_join(dat.lh, chromocounts)
dat.lh$mark <- sapply(dat.lh$experi, function(x) strsplit(x, "-")[[1]][[3]])

dat.lh <- dat.lh %>%
  dplyr::select(c(samp, chromocounts, spikeincounts, TA.frac))

dats.meta.annot.out2 <- lapply(dats.meta.annot.out, function(jdat){
  left_join(jdat, dat.lh, by = c("cell" = "samp"))
})


# Check bad cells from Buys plots -----------------------------------------

# check I filtered by number of cuts?


mats.bin <- lapply(mats, BinarizeMatrix)


# jmark <- "H3K4me3"
# inf.check <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_dropbox/from_slack/K562_3-134854697-138984717_", jmark, "_order.csv"))
# # inf.check <- file.path(hubprefix, "jyeung/data/scChiC/from_dropbox/from_slack/K562_3-134854697-138984717_H3K4me3_order.csv")
# # inf.check <- file.path(hubprefix, "jyeung/data/scChiC/from_dropbox/from_slack/K562_3-134854697-138984717_H3K27me3_order.csv")
# # inf.check <- file.path(hubprefix, "jyeung/data/scChiC/from_dropbox/from_slack/K562_3-134854697-138984717_H3K9me3_order.csv")
# dat.check <- fread(inf.check, skip = 1, header = TRUE)
# 
# rnames <- dat.check$V1
# dat.check <- as.data.frame(subset(dat.check, select = -V1))
# mat.check <- as.matrix(dat.check)
# rownames(mat.check) <- rnames
# 
# library(gplots)
# library(heatmap3)
# 
# # heatmap(t(mat.check))
# heatmap(mat.check)
# 
# mat.check.bin <- BinarizeMatrix(mat.check)
# 
# plot(mat.check.bin[1, ])
# plot(mat.check.bin[2, ])
# 
# nonzeros.frac <- apply(mat.check.bin, 1, function(jrow) nnzero(jrow) / length(jrow))
# 
# plot(density(nonzeros.frac))
# 
# jthres <- 1.5
# jsd <- mad(nonzeros.frac)
# jmean <- median(nonzeros.frac) 
# jmean + jsd * jthres
# 
# plot(density(nonzeros.frac))
# abline(v = jmean)
# abline(v = jmean + jsd *jthres, lty = 2)
# 
# bad.cells <- names(nonzeros.frac)[which(nonzeros.frac > jmean + jsd * jthres)]

nonzeros.global <- lapply(mats.bin, function(jmat){
  apply(jmat, 2, function(jcol) nnzero(jcol) / length(jcol))
})

bad.cells.lst <- lapply(jmarks, function(jmark){
  nonzeros.frac.tmp <- nonzeros.global[[jmark]]
  jsd.global <- mad(nonzeros.frac.tmp)
  jmean.global <- median(nonzeros.frac.tmp)
  jthres.frac <- jmean.global + jsd.global * jthres
  plot(density(nonzeros.frac.tmp), main = jmark, xlab = "Fraction of nonzero cuts globally")
  abline(v = jmean.global)
  abline(v = jthres.frac, lty = 2)
  bad.cells <- names(nonzeros.frac.tmp)[which(nonzeros.frac.tmp > jthres.frac)]
})

lapply(bad.cells.lst, length)

# Refilter ----------------------------------------------------------------


dats.meta.annot.out3 <- lapply(jmarks, function(jmark){
  bad.cells.tmp <- bad.cells.lst[[jmark]]
  jdat <- dats.meta.annot.out2[[jmark]] %>%
    ungroup() %>%
    filter(!cell %in% bad.cells.tmp)
})
lapply(dats.meta.annot.out2, dim)
lapply(dats.meta.annot.out3, dim)

dev.off()


# jdat <- data.frame(frac.nonzeros = nnzeros.global[[jmark]], cell = names(nnzeros.global[[jmark]]), stringsAsFactors = FALSE) %>%
#   ungroup() %>%
#   mutate(is.bad = cell %in% bad.cells)
# 
# ggplot(jdat, aes(x = frac.nonzeros, fill = is.bad)) +
#   geom_density(alpha = 0.25) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# range(subset(jdat, cell %in% bad.cells)$frac.nonzeros)
# 


# Check globally  ---------------------------------------------------------



# Write output ------------------------------------------------------------


for (jmark in jmarks){
  print(jmark)
  fname <- paste0("K562_cleaned_by_G1filt.", jmark, ".", jsuffix, ".txt")
  outf <- file.path(outdir, fname)
  # assertthat::assert_that(!file.exists(outf))
  fwrite(dats.meta.annot.out3[[jmark]], file = outf, sep = "\t")
}
