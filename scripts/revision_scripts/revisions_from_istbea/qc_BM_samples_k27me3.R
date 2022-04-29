t# Jake Yeung
# Date of Creation: 2022-01-25
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/qc_BM_samples.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

library(Matrix)
library(gtools)

frip.cutoff <- 0.8
ta.cutoff <- 0.5
counts.cutoff <- 3000

jmark <- "k27me3"
jmarkold <- "H3K27me3"
jexperi <- "BM"
inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/counts_tables_50000"
inmain.rz <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/LH_RZ_counts"
inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/LDA_outputs_first_submission/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_", jmarkold, ".cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_", jmarkold, ".cleaned.varfilt_2.K-30.Robj")

dname <- paste(jexperi, jmark, sep = "_")

outdir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA/", dname)
dir.create(outdir)

outf.allmerged <- file.path(outdir, paste0("count_mat_merged_with_old.", jmark, ".", Sys.Date(), ".rds"))
# outf.allmerged.fewer <- file.path(outdir, paste0("count_mat_merged_with_old_fewer2.", jmark, ".", Sys.Date(), ".rds"))
outf.newonly <- file.path(outdir, paste0("count_mat_new_only.", jmark, ".", Sys.Date(), ".rds"))
outf.qcmeta <- file.path(outdir, paste0("qc_metadata_new_only.", jmark, ".", Sys.Date(), ".txt"))
outpdf <- file.path(outdir, paste0("plots_", dname, ".", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)






# load count tables
indir <- file.path(inmain, dname)
indir.rz <- file.path(inmain.rz, dname)

# get inf mats
infs.fnames <- list.files(path = indir, pattern = "*.csv", full.names = FALSE)
infs.full <- list.files(path = indir, pattern = "*.csv", full.names = TRUE)
assertthat::assert_that(length(infs.fnames) > 0)

# get inf mats RZ LH 
infs.rz.fnames <- list.files(path = indir.rz, pattern = "*.csv", full.names = FALSE)
infs.rz.full <- list.files(path = indir.rz, pattern = "*.csv", full.names = TRUE)

names(infs.full) <- sapply(infs.fnames, function(f) strsplit(f, "\\.")[[1]][[1]])
names(infs.rz.full) <- sapply(infs.rz.fnames, function(f) strsplit(f, "\\.")[[1]][[1]])

mats.lst <- lapply(infs.full, function(jinf){
  ReadMatSlideWinFormat(jinf)
})

all.rnames <- unique(gtools::mixedsort(unlist(lapply(mats.lst, function(jmat) rownames(jmat)))))
mat.merged <- cbind.fill.lst(mats.lst, all.rnames = all.rnames, fill = 0)

dat.rz <- lapply(infs.rz.full, function(jinf){
  ReadLH.SummarizeTA(jinf)
}) %>%
  bind_rows()

dat.rz <- dat.rz %>%
  rowwise() %>%
  mutate(plate = ClipLast(samp, jsep = "_"))

# Plot total counts -------------------------------------------------------

ggplot(dat.rz, aes(x = log10(total.count), y = TA.frac)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz, aes(x = log10(total.count), y = TA.frac)) + 
  geom_point() + 
  facet_wrap(~plate) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot frips --------------------------------------------------------------

bin.avgs <- log10(rowMeans(as.matrix(mat.merged)) * 1000  + 1)

# define peaks arbitrarily
bin.cutoff <- mean(bin.avgs)
peaks <- names(bin.avgs)[bin.avgs > bin.cutoff]

plot(density(bin.avgs))
abline(v = bin.cutoff)

# Calculate fraction of cells in or out of peaks  -------------------------

cell.sums <- Matrix::colSums(mat.merged)
cell.sums.atpeaks <- Matrix::colSums(mat.merged[peaks, ])

cell.data <- data.frame(cell = colnames(mat.merged), total.count.from.mat = cell.sums, count.at.peak = cell.sums.atpeaks, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(frac.count.in.peak = count.at.peak / total.count.from.mat)



dat.rz.annot <- left_join(dat.rz, cell.data, by = c("samp" = "cell"))

ggplot(dat.rz.annot, aes(x = frac.count.in.peak)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  geom_vline(xintercept = frip.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = frac.count.in.peak)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~plate) + 
  theme_bw() +
  geom_vline(xintercept = frip.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = total.count)) + 
  geom_density() + 
  theme_bw() + 
  scale_x_log10() + 
  geom_vline(xintercept = counts.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = total.count)) + 
  geom_density() + 
  theme_bw() + 
  scale_x_log10() + 
  facet_wrap(~plate) + 
  geom_vline(xintercept = counts.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = TA.frac)) + 
  geom_density() + 
  theme_bw() + 
  geom_vline(xintercept = ta.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = TA.frac)) + 
  geom_density() + 
  theme_bw() + 
  facet_wrap(~plate) + 
  geom_vline(xintercept = ta.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# select good cells 

dat.rz.annot <- dat.rz.annot %>%
  rowwise() %>%
  mutate(is.good = frac.count.in.peak >= frip.cutoff & total.count >= counts.cutoff & TA.frac >= ta.cutoff)

dat.rz.annot.filt <- subset(dat.rz.annot, is.good)

ggplot(dat.rz.annot, aes(x = total.count, y = TA.frac, color = is.good)) + 
  geom_point() + 
  theme_bw() + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = frac.count.in.peak, fill = is.good)) + 
  geom_density(alpha = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = total.count, y = TA.frac, color = is.good)) + 
  geom_point() + 
  facet_wrap(~plate) + 
  theme_bw() + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Load H3K27me3 counts ----------------------------------------------------

load(inf.lda, v=T)
count.mat.orig <- count.mat

# get common rows
rows.common <- intersect(rownames(count.mat.orig), rownames(mat.merged))
assertthat::assert_that(length(rows.common) > 0)


# Merge together ----------------------------------------------------------

cnames.orig <- colnames(count.mat.orig)
cnames.plates <- unique(sapply(cnames.orig, function(x) ClipLast(x, jsep = "_")))
# cnames.orig.fewer <- cnames.orig[grepl("rep2|rep3", cnames.orig)]
# assertthat::assert_that(length(cnames.orig.fewer) > 0)


count.mat.merged.with.old <- cbind(count.mat.orig[rows.common, ], mat.merged[rows.common, ])
# mat.merged <- mat.merged[rows.common, ]
# count.mat.merged.with.old <- cbind(count.mat.orig[rows.common, cnames.orig.fewer], mat.merged[rows.common, ])

# Write outputs -----------------------------------------------------------

saveRDS(count.mat.merged.with.old, outf.allmerged)
# saveRDS(count.mat.merged.with.old.fewer, outf.allmerged.fewer)
saveRDS(mat.merged[rows.common, ], file = outf.newonly)
fwrite(dat.rz.annot, file = outf.qcmeta, sep = "\t")

dev.off()

