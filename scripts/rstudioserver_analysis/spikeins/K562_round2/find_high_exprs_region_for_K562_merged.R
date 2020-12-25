# Jake Yeung
# Date of Creation: 2020-11-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/find_high_exprs_region_for_K562_merged.R
# # calculate sensitivity of the assay 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

# Get raw -----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmark <- "H3K4me3"
print(jmark)


inf.raw <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters/K562_AllMerged_", jmark, ".merged.sorted.tagged.countTable.binsize_20000.csv")
dat.raw <- ReadMatSlideWinFormat(inf.raw)

# good cells ----------------------------------------------------------------

# jmark <- "H3K27me3"
indir.good <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2.G1filt")
outdir <- indir.good
# infs.good <- list.files()
fname <- paste0("K562_count_tables_50000.", jmark, ".G1filt.rds")
# outfname <- paste0("K562_merged_good_cells_filt_20000.", jmark, ".G1filt.bed")
# outpath <- file.path(outdir, outfname)
# assertthat::assert_that(!file.exists(outpath))

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

# search chr3:136,463,486-137,065,453
jbin <- grep("chr3:13618", names(jcounts.vec), value = TRUE)[[1]]
jcounts.vec[[jbin]]

plot(density(jcounts.vec), main = "Distribution of # cuts in 1075 K562 cells for K4me3. 20kb bins", log = "x")


plot(density(jcounts.vec))
