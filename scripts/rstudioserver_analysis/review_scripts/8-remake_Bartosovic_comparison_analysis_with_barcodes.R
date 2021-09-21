# Jake Yeung
# Date of Creation: 2021-06-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/8-remake_Bartosovic_comparison_analysis_with_barcodes.R
# Copied from https://github.com/Castelo-Branco-lab/scCut-Tag_2020/blob/949c19a44a6e98444d40e5b8cd253f17d2e54aad/notebooks/revision_unsorted_code.Rmd#L42-L193 


rm(list=ls())

# library(tidyr)
# library(data.table)
# library(Matrix)


# UMAP
library(dplyr)
library(data.table)
library(ggplot2)

library(Seurat)
library(ggplot2)
library(ggthemes)
library(viridis)
# library(Signac)
library(reshape2)
library(scales)
library(reticulate)
library(ggrepel)

params <- list()
params$out_prefix <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic/scCT"
assertthat::assert_that(dir.exists(params$out_prefix))


######### FriP analysis ploting

prefix <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic"

# prefix   <- paste0(params$out_prefix,'other_datasets/frip_analysis/')
files    <- c('scCT/H3K27me3_N1/','scCT/H3K27me3_N2/','scCT/H3K27me3_N3/','scCT/H3K27me3_N4/','kaya_okur/K562_H3K4me2_iCell8/','kaya_okur/K562_H3K27me3_iCell8/','kaya_okur/H1_H3K27me3_iCell8/','scCT/H3K27me3_cell_lines_1/','scCT/H3K27me3_cell_lines_2/')
# files    <- c('scCT/H3K27me3_N1/','scCT/H3K27me3_N2/','scCT/H3K27me3_N3/','scCT/H3K27me3_N4/','kaya_okur/K562_H3K4me2_iCell8/','kaya_okur/K562_H3K27me3_iCell8/','kaya_okur/H1_H3K27me3_iCell8/','scCT/H3K27me3_cell_lines_1/','scCT/H3K27me3_cell_lines_2/')
group    <- c(rep("scCT_brain",4),rep("Kaya-Okur",3),rep('scCT_cell_lines',2))
sample   <- gsub("scCT|kaya_okur|/","",files)
antibody <- c(rep('H3K27me3',4),'H3K4me2',rep('H3K27me3',4))
all.sufix   <- 'all_fragments.txt'
peaks.sufix <- 'peak_fragments.txt'

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load barcodes  ----------------------------------------------------------

inf.primary.bcs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic/H3K27me3_cluster_barcode_table.csv"
inf.cell_lines.bcs <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic/H3K27me3_cell_lines_cluster_barcode_table.csv"

dat.primary.bcs <- fread(inf.primary.bcs, header = TRUE, col.names = c("cluster", "cell"))
dat.cell_lines.bcs <- fread(inf.cell_lines.bcs, header = TRUE, col.names = c("cluster", "cell"))


# Load cells --------------------------------------------------------------

# inf.barcodes <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds"
# dat.barcodes <- readRDS(inf.barcodes)
# dat.barcodes.bart <- subset(dat.barcodes, grepl("^Bart", jset))

dat.barcodes.bart <- bind_rows(dat.primary.bcs, dat.cell_lines.bcs)

# hubprefix <- "jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic/scCT/H3K27me3_cell_lines_1"
inf.allfrags.bart <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic/scCT/H3K27me3_cell_lines_1/all_fragments.txt")
inf.peakfrags.bart <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic/scCT/H3K27me3_cell_lines_1/peak_fragments.txt")

dat.allfrags.bart <- fread(file.path(inf.allfrags.bart))
dat.peakfrags.bart <- fread(file.path(inf.peakfrags.bart))

dbase <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Jake_request_from_MarekBartosovic")
dvec <- list("scCT/H3K27me3_N1", "scCT/H3K27me3_N2", "scCT/H3K27me3_N3", "scCT/H3K27me3_N4", "scCT/H3K27me3_cell_lines_1", "scCT/H3K27me3_cell_lines_2", "kaya_okur/H1_H3K27me3_iCell8", "kaya_okur/K562_H3K27me3_iCell8", "grosselin/SRR7536861", "grosselin/SRR7536860", "grosselin/SRR7536862")
bnames <- unlist(lapply(dvec, function(x) dirname(x)))
dnames <- unlist(lapply(dvec, function(x) x))
names(dnames) <- dnames
names(dvec) <- names(dnames)
names(bnames) <- names(dnames)

infs.allfrags <- lapply(dnames, function(dname){
  d <- dvec[[dname]]
  indir <- file.path(dbase, d)
  fname <- list.files(indir, pattern = "*all*")
  inf.tmp <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

infs.peakfrags <- lapply(dnames, function(dname){
  d <- dvec[[dname]]
  indir <- file.path(dbase, d)
  fname <- list.files(indir, pattern = "*peak*")
  inf.tmp <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.allfrags <- lapply(dnames, function(dname){
  dat <- fread(infs.allfrags[[dname]], col.names = c("allcounts", "cell"))
  # dat$type <- "allfrags"
  dat$jset <- dname
  dat$paper <- bnames[[dname]]
  return(dat)
}) %>%
  bind_rows()

dat.peakfrags <- lapply(dnames, function(dname){
  dat <- fread(infs.peakfrags[[dname]], col.names = c("peakcounts", "cell"))
  # dat$type <- "allfrags"
  dat$jset <- dname
  dat$paper <- bnames[[dname]]
  return(dat)
}) %>%
  bind_rows()

dat.frags.merged <- left_join(dat.peakfrags, dat.allfrags)

# remove some barcodes? 
barcodes.keep <- unique(dat.barcodes.bart$cell)
barcodes.remove <- subset(dat.frags.merged, paper == "scCT" & !cell %in% barcodes.keep)$cell

dat.frags.merged.filt <- subset(dat.frags.merged, !cell %in% barcodes.remove)

ggplot(dat.frags.merged.filt, aes(x = jset, y = peakcounts / allcounts)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data"
outname <- paste0("frip_peakcounts_allcounts.all_scCT.", Sys.Date(), ".rds")
outrds <- file.path(outdir, outname)
saveRDS(dat.frags.merged.filt, file = outrds)

# inf.test <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1_notrim_bowtie2/hiddendomains_output/SRR12638101_1.sorted.1000.cutoff/SRR12638101_1.sorted.1000.cutoff_treatment_bins.txt")
# dat.test <- fread(inf.test)
