# Jake Yeung
# Date of Creation: 2022-01-09
# File: ~/projects/scchic/scripts/macbook_analysis_2021/new_experiments/qc_BM_samples.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


jmarks <- c("k27me3", "k9me3", "k4me1-k9me3"); names(jmarks) <- jmarks
jmark <- "k27me3"
jexperi <- "BM"
inmain <- "/Users/yeung/data/scchic/from_cluster_2021/new_experiments/count_tables/counts_tables_50000"
inmain.rz <- "/Users/yeung/data/scchic/from_cluster_2021/new_experiments/count_tables/LH_RZ_counts"

dname <- paste(jexperi, jmark, sep = "_")

outdir <- paste0("/Users/yeung/data/scchic/from_cluster_2021/new_experiments/filtered_count_tables_for_LDA/", dname)
dir.create(outdir)

outf.allmerged <- file.path(outdir, paste0("count_mat_merged_with_old.", jmark, ".", Sys.Date(), ".rds"))
outf.newonly <- file.path(outdir, paste0("count_mat_new_only.", jmark, ".", Sys.Date(), ".rds"))
outf.qcmeta <- file.path(outdir, paste0("qc_metadata_new_only.", jmark, ".", Sys.Date(), ".txt"))
outpdf <- file.path(outdir, paste0("plots_", dname, ".", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# Do stuff ----------------------------------------------------------------




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

library(Matrix)
library(gtools)
mats.lst <- lapply(infs.full, function(jinf){
  ReadMatSlideWinFormat(jinf)
})

all.rnames <- gtools::mixedsort(unlist(lapply(mats.lst, function(jmat) rownames(jmat))))
mat.merged <- cbind.fill.lst(mats.lst, all.rnames = all.rnames, fill = 0)

dat.rz <- lapply(infs.rz.full, function(jinf){
  ReadLH.SummarizeTA(jinf)
}) %>%
  bind_rows()


# Plot total counts -------------------------------------------------------

ggplot(dat.rz, aes(x = log10(total.count), y = TA.frac)) + 
  geom_point() + 
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

frip.cutoff <- 0.8
ta.cutoff <- 0.5
counts.cutoff <- 3000


dat.rz.annot <- left_join(dat.rz, cell.data, by = c("samp" = "cell"))

ggplot(dat.rz.annot, aes(x = frac.count.in.peak)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  geom_vline(xintercept = frip.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = total.count)) + 
  geom_density() + 
  theme_bw() + 
  scale_x_log10() + 
  geom_vline(xintercept = counts.cutoff, linetype = "dotted") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.rz.annot, aes(x = TA.frac)) + 
  geom_density() + 
  theme_bw() + 
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


# Load H3K27me3 counts ----------------------------------------------------

inf.lda <- "/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"
load(inf.lda, v=T)
count.mat.orig <- count.mat

# get common rows
rows.common <- intersect(rownames(count.mat.orig), rownames(mat.merged))
assertthat::assert_that(length(rows.common) > 0)


# Merge together ----------------------------------------------------------

count.mat.merged.with.old <- cbind(count.mat.orig[rows.common, ], mat.merged[rows.common, ])


# Write outputs -----------------------------------------------------------


saveRDS(count.mat.merged.with.old, outf.allmerged)
saveRDS(mat.merged[rows.common, ], file = outf.newonly)
fwrite(dat.rz.annot, file = outf.qcmeta, sep = "\t")

dev.off()

