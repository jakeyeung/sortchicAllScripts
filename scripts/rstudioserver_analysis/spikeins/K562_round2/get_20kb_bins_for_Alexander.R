# Jake Yeung
# Date of Creation: 2020-10-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/get_20kb_bins_for_Alexander.R
# Merge 20kb bins to ocmpare with cime

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

# Get raw -----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  print(jmark)
  
  
  inf.raw <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters/K562_AllMerged_", jmark, ".merged.sorted.tagged.countTable.binsize_20000.csv")
  dat.raw <- ReadMatSlideWinFormat(inf.raw)
  
  # good cells ----------------------------------------------------------------
  
  # jmark <- "H3K27me3"
  indir.good <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2.G1filt")
  outdir <- indir.good
  # infs.good <- list.files()
  fname <- paste0("K562_count_tables_50000.", jmark, ".G1filt.rds")
  outfname <- paste0("K562_merged_good_cells_filt_20000.", jmark, ".G1filt.bed")
  outpath <- file.path(outdir, outfname)
  assertthat::assert_that(!file.exists(outpath))
  
  inf.good <- file.path(indir.good, fname)
  assertthat::assert_that(file.exists(inf.good))
  dat.good <- readRDS(inf.good)
  
  
  # Filter good cells, merge, write to file  --------------------------------
  
  good.cells <- colnames(dat.good)
  
  cnames.keep <- which(colnames(dat.raw) %in% good.cells)
  assertthat::assert_that(length(cnames.keep) > 0)
  dat.raw.filt <- dat.raw[, cnames.keep]
  
  coords <- rownames(dat.raw.filt)
  jcounts.vec <- rowSums(dat.raw.filt)
  outbed <- data.frame(chromo = sapply(coords, GetChromo), 
                       jstart = sapply(coords, GetStart),
                       jend = sapply(coords, GetEnd),
                       jcounts = jcounts.vec,
                       stringsAsFactors = FALSE)
  
  fwrite(outbed, file = outpath, sep = "\t", col.names = FALSE, row.names = FALSE)
}



