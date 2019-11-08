# Jake Yeung
# Date of Creation: 2019-10-22
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/1-get_good_cells_OUD3909_OUD3910.R

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)


# Functions ---------------------------------------------------------------

ReadRZ <- function(inf, remove.nones = FALSE){
  dat.tmp <- as.data.frame(fread(inf, sep = ",", header = TRUE))
  colnames(dat.tmp)[[1]] <- dat.tmp[1, 1]
  dat.tmp <- dat.tmp[-1, ]
  
  if (remove.nones){
    dat.tmp <- subset(dat.tmp, recognizedSequence != "None")
  }
  
  dat.mat <- as.matrix(dat.tmp[, -1])
  dat.mat.sum <- colSums(dat.mat, na.rm = TRUE)
  dat.mat.ta <- subset(dat.tmp, recognizedSequence == "TA", select = -recognizedSequence)  # I think Buys calls it AT for some reason here? Need to check actually, might be an old bam file
  dat.merged <- rbind(dat.mat.ta, dat.mat.sum)
  # add first column
  dat.merged.annot <- cbind(data.frame(V1 = c("TA_start", "total")), dat.merged)  # call it V1 and rownames same as Buys TA_obs_per_cell.csv
  return(dat.merged.annot)
}


# Outputs -----------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells"
outdir <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells_part2"
dir.create(outdir)

# Define files ------------------------------------------------------------

# first the RZ
ouds <- c("oud3909", "oud3910")
names(ouds) <- ouds

jsuffix <- ".Ensembl98"
# jsuffix <- ""

if (jsuffix == ""){
  jgenome <- "Ensembl93"
} else {
  jgenome <- "Ensembl98"
}

indir.rz.lst <- lapply(ouds, function(oud) paste0("/Users/yeung/data/scchic/from_cluster/", oud, "/", oud, "_QC_plots", jsuffix, "/", oud, "_HD_0/RZcounts"))
infs.rz <- unlist(lapply(indir.rz.lst, function(indir.rz) list.files(indir.rz, pattern = "PZ-ChIC-ZFWKM.*RZ_counts.csv", full.names = TRUE)))
assertthat::assert_that(all(sapply(indir.rz.lst, function(indir.rz) dir.exists(indir.rz))))
assertthat::assert_that(all(file.exists(infs.rz)))

# count tables
indir.count.lst <- lapply(ouds, function(oud) paste0("/Users/yeung/data/scchic/from_cluster/", oud, "/", oud, "_QC_plots", jsuffix, "/", oud, "_HD_0/countTables"))
infs.count <- unlist(lapply(indir.count.lst, function(indir.count) list.files(indir.count, pattern = "PZ-ChIC-ZFWKM.*countTable.rds", full.names = TRUE)))
assertthat::assert_that(all(sapply(indir.count.lst, function(count.rz) dir.exists(count.rz))))
assertthat::assert_that(all(file.exists(infs.count)))

jnames <- lapply(infs.count, function(inf){
  bname <- basename(inf)
  bname <- sapply(bname, function(x) strsplit(x, "\\.")[[1]][[1]])
})
names(infs.count) <- jnames

# Load data ---------------------------------------------------------------


dats <- lapply(infs.rz, function(inf){
  dat <- ReadRZ(inf, remove.nones = TRUE)
  # make long
  dat.TA <- tidyr::gather(dat %>% filter(V1 == "TA_start"), key = "samp", value = "count", -V1) %>%
    dplyr::rename(TA.count = count) %>%
    dplyr::select(-V1)
  dat.Total <- tidyr::gather(dat %>% filter(V1 == "total"), key = "samp", value = "count", -V1) %>%
    dplyr::rename( total.count = count) %>%
    dplyr::select(-V1)
  # merge the two
  dat.long <- left_join(dat.TA, dat.Total, by = c("samp"))
}) %>%
  bind_rows() %>%
  mutate(TA.frac = TA.count / total.count)

# Get Index and Experi
dats$experi <- sapply(dats$samp, function(x) strsplit(x, "_")[[1]][[1]])
dats$indx <- sapply(dats$experi, function(x){
  jsplit <- strsplit(x, "-")[[1]]
  return(jsplit[length(jsplit) - 3])
})
dats$is.stem.enriched <- sapply(dats$experi, function(x){
  return(grepl("stem-cell", x))
})
dats$cell <- sapply(dats$samp, function(x){
  return(paste0("cell", strsplit(x, "_")[[1]][[2]]))
})

# check samp names are unique
assertthat::assert_that(length(which(duplicated(subset(dats)$samp))) == 0)

ggplot(dats, aes(x = log10(total.count), y = TA.frac)) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi) + 
  ggtitle(jgenome)


empty.wells <- GetEmptyWells()

dats$is.empty <- sapply(dats$cell, function(x) x %in% empty.wells)

# set a cutoff: some bug in TA vs total count? need to ask Buys
jcutoff <- 2.5
ta.cutoff <- 0.5

dats <- dats %>%
  rowwise() %>%
  mutate(good.cell = (TA.frac >= ta.cutoff & log10(total.count) >= jcutoff))

m <- ggplot(dats, aes(x = log10(total.count), y = TA.frac, color = good.cell, shape = is.empty)) + geom_point() +
  theme_bw(6) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi) + 
  geom_vline(xintercept = jcutoff) + 
  geom_hline(yintercept = ta.cutoff) + 
  ggtitle(jgenome)

pdf(paste0("/Users/yeung/data/scchic/quality_control_ZF/mats_filt_by.good_cells.sizecutoff_", jcutoff, ".TAcutoff_", ta.cutoff, ".pdf"))
print(m)
dev.off()


cells.keep <- subset(dats, good.cell)$samp


# Remove bad plates -------------------------------------------------------

dats.summary <- dats %>%
  group_by(experi) %>%
  summarise(ncell = length(which(good.cell == TRUE)))
plates.keep <- subset(dats.summary, ncell > 0)$experi

infs.count.filt <- infs.count[which(names(infs.count) %in% plates.keep)]

# Write good cells --------------------------------------------------------

# keep all bins and filter them out later 
mats <- lapply(infs.count.filt, function(inf){
  mat.tmp <- readRDS(inf)
  # rows.keep <- which(rownames(mat.tmp) %in% jterms)
  # return(mat.tmp[rows.keep, ])
})

lapply(mats, function(x) nnzero(x) / length(x))
lapply(mats, dim)

# # get common rows. By mark?
# jmarks <- unique(sapply(names(mats), function(x) strsplit(x, "-")[[1]][[4]], USE.NAMES = FALSE))
# names(jmarks) <- jmarks

# mats.filt <- mats[grepl(jmark, names(infs.count))]
# mats.bymark <- lapply(jmarks, function(jmark){
#   return(mats.filt)
# })

rnames <- lapply(mats, rownames)
rnames.common <- Reduce(intersect, rnames)

print(length(rnames.common))

# rnames.bymark <- lapply(mats.bymark, function(mats.lst) lapply(mats.lst, rownames))
# rnames.common.bymark <- lapply(rnames.bymark, function(rnames) Reduce(intersect, rnames))
# lapply(rnames.common.bymark, function(x) length(x))
# 
# do one more filter and write to output

mats.filt <- lapply(mats, function(x){
  rows.keep <- which(rownames(x) %in% rnames.common)
  cols.keep <- which(colnames(x) %in% cells.keep)
  return(x[rows.keep, cols.keep])
})
lapply(mats.filt, dim)


# Write all  --------------------------------------------------------------

# write as list, keep things simpler later when we filter for good bins

# write outputs

outf <- paste0("/Users/yeung/data/scchic/quality_control_ZF/mats_filt_by.good_cells.sizecutoff_", jcutoff, ".TAcutoff_", ta.cutoff, ".rds")
saveRDS(mats.filt, file = outf)


