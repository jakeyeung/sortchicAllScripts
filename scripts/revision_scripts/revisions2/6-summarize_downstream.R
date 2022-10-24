# Jake Yeung
# Date of Creation: 2022-07-25
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/6-summarize_downstream.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks


# Cells keep --------------------------------------------------------------

inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})


# Load RZs ----------------------------------------------------------------


indir.rz <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/LH_RZ_counts_longer"

jmarkold2new <- hash::hash(jmarksold, names(jmarksold))
dat.rz.lst <- lapply(jmarksold, function(jmark){
  print(jmark)
  inf <- file.path(indir.rz, paste0("BM_allmerged_", jmark, ".LH_counts.csv"))
  jmarknew <- jmarkold2new[[jmark]]
  cells.keep <- dat.meta.lst[[jmarknew]]$cell
  dat.rz <- scchicFuncs::ReadLH.SummarizeTA(inf) %>% 
    mutate(mark = jmark) %>%
    filter(samp %in% cells.keep)
})

# try loading 50kb genome wide bins instead
indir.bins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_bins/counts_tables_50000"
total.cuts.lst <- lapply(jmarksold, function(jmark){
  print(jmark)
  indir.tmp <- file.path(indir.bins, paste0("BM_", jmark))
  inf.tmp <- file.path(indir.tmp, paste0("BM_allmerged_", jmark, ".countTable.binsize_50000.csv"))
  mat.bins <- ReadMatSlideWinFormat(inf.tmp)
  cell.sums <- colSums(mat.bins)
  dat.cell.sums <- data.frame(samp = names(cell.sums), cut.total.allbins = cell.sums, stringsAsFactors = FALSE)
})


inf.mat.k9 <- file.path(indir.bins, paste0("BM_H3K9me3"), paste0("BM_allmerged_H3K9me3.countTable.binsize_50000.csv"))
count.mat.k9 <- ReadMatSlideWinFormat(inf.mat.k9)
pbulk.k9 <- rowSums(count.mat.k9)
pbulk.k9.norm <- pbulk.k9 / sum(pbulk.k9)
plot(density(pbulk.k9.norm))
bg.log10 <- 0.000015

plot(density(log10(pbulk.k9.norm + 1)))
abline(v = bg.log10)

# get dynamic k9 bins 
inf.k9.bins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/dynamic_bins/dynamic_bins_50kb.k9me3.2022-07-24.txt"
dat.k9.bins <- fread(inf.k9.bins)$V4

dat.mat.k9.annot <- data.frame(signal.linear = pbulk.k9.norm, bin = names(pbulk.k9.norm)) %>%
  rowwise() %>%
  mutate(is.dyn = bin %in% dat.k9.bins)

ggplot(dat.mat.k9.annot, aes(x = log10(signal.linear + 1), fill = is.dyn)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept = bg.log10) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add summary for different bins  -----------------------------------------

inf.dynbins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/counts_different_bins/summary_bins_objs_allmarks_dynbinslst.2022-07-25.RData"
load(inf.dynbins, v=T)

inf.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/counts_different_bins/summary_bins_objs_TSS.2022-07-25.RData"
load(inf.tss, v=T)

head(summary.dynbins.bybin.lst$k4me3$k4me3)
head(summary.tss.lst$k4me3)

# edit colnames for easier joining later 
summary.dynbins.renamed.lst <- lapply(jmarks, function(jmark1){
  summary.dynbins.bybin.lst <- lapply(jmarks, function(jmark2){
    dat.tmp <- summary.dynbins.bybin.lst[[jmark1]][[jmark2]]
    # colnames "samp"       "cuts.total" "jtype"      "nbins"
    suffix <- jmark2
    cnames.old <- colnames(dat.tmp)[2:ncol(dat.tmp)]
    cnames.new <- paste(cnames.old, suffix, sep = ".")
    colnames(dat.tmp) <- c("samp", cnames.new)
    return(dat.tmp)
  })
})

# merge them all together into one: left joins

summary.tss.renamed.lst <- lapply(jmarks, function(jmark){
  dat.tss <- summary.tss.lst[[jmark]] %>%
    dplyr::rename(cuts.total.TSS = cuts.total, jtype.TSS = jtype, nbins.TSS = nbins)
  return(dat.tss)
})

dat.summary.all <- lapply(jmarks, function(jmark){
  dat.joined <- Reduce(f = left_join,x = summary.dynbins.renamed.lst[[jmark]], init = summary.tss.renamed.lst[[jmark]])
})


# Combine with RZ  ---------------------------------------------------------

suffixs <- c("TSS", "k4me1", "k4me3", "k27me3", "k9me3"); names(suffixs) <- suffixs
# normalize by total and by number of bins
dat.summary.all.with.total <- lapply(jmarks, function(jmark){
  print(jmark)
  jdat <- left_join(dat.rz.lst[[jmark]], dat.summary.all[[jmark]]) %>%
    rowwise()
  for (suffix in suffixs){
    print(suffix)
    totalname <- "total.count"
    newname <- paste0("cuts.", suffix, ".norm")
    newname.bin <- paste0("cuts.", suffix, ".norm.bybins")
    cutname <- paste0("cuts.total.", suffix)
    bname <- paste0("nbins.", suffix)
    jdat[[newname]] <- jdat[[cutname]] / jdat[[totalname]]
    jdat[[newname.bin]] <- jdat[[newname]] / jdat[[bname]]
    # cnames.keep <- c("samp", newname, newname.bin, "total.count")
    # jdat.filt <- jdat[, cnames.keep]
  }
  # filter colnames
  cnames.keep1 <- paste0("cuts.", suffixs, ".norm")
  cnames.keep2 <- paste0("cuts.", suffixs, ".norm.bybins")
  cnames.keep <- c("samp", cnames.keep1, cnames.keep2, "total.count")
  jdat.filt <- jdat[, cnames.keep]
  # make long
  jdat.long <- melt(jdat.filt, id.vars = c("samp", "total.count"), variable.name = "jtype", value.name = "normval") %>%
    rowwise() %>%
    mutate(suffix = ifelse(endsWith(x = as.character(jtype), ".norm"), "bytotal", "bybins"))
  jdat.long.lst <- split(jdat.long, f = jdat.long$suffix)
  return(jdat.long.lst)
})

# # checks
# jmarkcheck <- "k9me3"
# jcheck <- dat.summary.all.with.total[[jmarkcheck]] %>%
#   mutate(check = cuts.total.TSS / total.count)
# 
# ggplot(jcheck, aes(x = check)) +
#   geom_histogram() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 

# Summarize fraction across different types of bins -----------------------

m.lst <- lapply(jmarks, function(jmarkcheck){
  m <- ggplot(dat.summary.all.with.total[[jmarkcheck]]$bybins, aes(x = jtype, y = normval)) + 
    geom_boxplot() + 
    theme_bw() + 
    ggtitle(jmarkcheck) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
})
print(m.lst)


# Check by TSS norm across marks ------------------------------------------

dat.summary.long <- lapply(jmarks, function(jmark){
  jdat <- dat.summary.all.with.total[[jmark]]$bybins %>%
  # jdat <- dat.summary.all.with.total[[jmark]]$bytotal %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()


m.tss <- ggplot(dat.summary.long, aes(x = mark, y = normval)) + 
  geom_boxplot() + 
  facet_wrap(~jtype) + 
  theme_bw() + 
  # geom_hline(yintercept = 10^bg.log10 - 1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.tss)

m.tss.log <- ggplot(dat.summary.long, aes(x = mark, y = log10(normval + 1))) + 
  geom_boxplot() + 
  facet_wrap(~jtype) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.tss.log)







# 
# indir.sum <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/counts_different_bins"
# 
# jmarktmp <- "H3K4me1"
# 
# summary.b
# inf.tmp <- file.path(indir.sum, paste0("summary_bins_objs_", jmarktmp, "_dynbinslst.2022-07-25.RData"))
# load(inf.tmp, v=T)
# head(summary.dynbins.lst$k4me1)
# 



