# Jake Yeung
# Date of Creation: 2020-10-18
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/filter_good_cells_for_bam_clustering.G1only.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)


# Load mats ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2.G1filt"

indir.chromo <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
indir.lh <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_TAfrac.NewFilters")


mats <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("K562_count_tables_50000.", jmark, ".G1filt.rds")
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

# Write output ------------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K562_clean_by_G1filt"
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  fname <- paste0("K562_cleaned_by_G1filt.", jmark, ".txt")
  outf <- file.path(outdir, fname)
  assertthat::assert_that(!file.exists(outf))
  fwrite(dats.meta.annot.out2[[jmark]], file = outf, sep = "\t")
}
