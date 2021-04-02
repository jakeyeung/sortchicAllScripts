# Jake Yeung
# Date of Creation: 2021-03-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/dendrograms/3-calculate_variance_on_pseudobulk.CPM.from_50kb_bins.no_heatmap.mean_deviation_from_HSPCs.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)

library(hash)

jymax <- 25000
jymaxnorm <- 1.25
jymaxfrac <- 0.5

make.plots <- FALSE
# make.plots <- TRUE

# merge.ctypes.by.lineage <- FALSE
# merge.ctypes.by.lineage <- FALSE

jscale <- 1000000L
# log10filt <- -2
log10filt <- -1

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps"
outpdf <- file.path(outdir, paste0("variance_by_clusters_by_HighBins.", Sys.Date(), ".", jscale, ".log10filt_", log10filt, ".noheatmap.bootstraps.deviations.pdf"))

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

# merge some celltypes
ctypes <- list("Eryths" = "Erythroid", 
               "Bcells" = "Lymphoid", 
               "NKs" = "Lymphoid", 
               "Granulocytes" = "Myeloid",
               "Basophils" = "Myeloid", 
               "pDCs" = "Lymphoid",
               "DCs" = "Myeloid", 
               "HSPCs" = "HSPCs",
               "Erythroid" = "Erythroid",
               "Lymphoid" = "Lymphoid",
               "Myeloid" = "Myeloid")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load metadata -----------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(file.path(indir.meta, fname)) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[cluster]])
}) 
 
# cluster to col
cluster2col <- hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
cluster2col[["Erythroid"]] <- "#0072B2"
cluster2col[["Lymphoid"]] <- "#56B4E9"
cluster2col[["Myeloid"]] <- "#D55E00"

cname2color <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
colpalette <- colorRampPalette(c("grey1", "grey35", "grey99"))(1024)

# Load mats ---------------------------------------------------------------


inf.lda.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  assertthat::assert_that(file.exists(inf.lda))
  return(inf.lda)
})

count.mat.bins.raw.lst <- lapply(inf.lda.lst, function(jinf){
  load(jinf, v=T)
  return(count.mat)
})

# 
# 
# inf.lda.bins.lst <- lapply(jmarks, function(jmark){
#   if (jmark != "H3K27me3"){
#     inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
#   } else {
#     inf.lda.tmp <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_", jmark, ".cleaned.varfilt_2.K-30.Robj"))
#   }
#   return(inf.lda.tmp)
# })
# 
# 
# count.mat.bins.lst <- lapply(inf.lda.bins.lst, function(inf.bins){
#   print(inf.bins)
#   inf.bins <- inf.bins
#   load(inf.bins, v=T)  # out.lda, count.mat
#   return(count.mat)
# })
# 


# Load high bins ----------------------------------------------------------

bsize <- 50000
dat.high.bins <- lapply(jmarks, function(jmark){
  indir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned")
  inf <- file.path(indir, paste0("High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt"))
  dat <- fread(inf)
})

high.bins.lst <- lapply(dat.high.bins, function(jdat){
  jdat$CoordOriginal
})

# check 
assertthat::assert_that(!all(unlist(lapply(high.bins.lst, function(x) any(grepl("\\+", x))))))


# Filter bins -------------------------------------------------------------

count.mat.peaks.raw.filt.lst <- lapply(jmarks, function(jmark){
  count.mat <- count.mat.bins.raw.lst[[jmark]]
  # count.mat <- count.mat.bins.lst[[jmark]]
  # rmeans <- log10(rowMeans(count.mat))
  bfilt.keep <- high.bins.lst[[jmark]]
  # bfilt.keep <- bins.high.lst[[jmark]]  # from File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/14-check_differential_bins_and_peaks.R
  bfilt.keep2 <- which(rownames(count.mat) %in% bfilt.keep)
  assertthat::assert_that(length(bfilt.keep) == length(bfilt.keep2))
  print(length(bfilt.keep))
  print(length(bfilt.keep2))
  # bfilt.keep2.check <- bfilt.keep[which(!bfilt.keep %in% rownames(count.mat))]
  count.mat[bfilt.keep, ]
})

lapply(count.mat.bins.raw.lst, dim)
lapply(count.mat.peaks.raw.filt.lst, dim)


# Filter by rmaens? -------------------------------------------------------

rmeans.lst <- lapply(count.mat.peaks.raw.filt.lst, function(x){
  log10(rowMeans(x))
})

for (jmark in jmarks){
  plot(density(rmeans.lst[[jmark]]), main = jmark)
}

count.mat.peaks.raw.filt.lst2 <- lapply(jmarks, function(jmark){
  jmat.tmp <- count.mat.peaks.raw.filt.lst[[jmark]]
  rmeans <- rmeans.lst[[jmark]]
  rows.keep <- rmeans > log10filt
  jmat.tmp2 <- jmat.tmp[rows.keep, ]
})

lapply(count.mat.peaks.raw.filt.lst, dim)
lapply(count.mat.peaks.raw.filt.lst2, dim)

count.mat.peaks.raw.filt.lst <- count.mat.peaks.raw.filt.lst2

# Get pbu,k ---------------------------------------------------------------


dat.raw.pbulk.lst <- lapply(jmarks, function(jmark){
  count.mat <- count.mat.peaks.raw.filt.lst[[jmark]]
  cnames.keep.lst <- split(x = dat.metas[[jmark]]$cell, f = dat.metas[[jmark]]$cluster)
  dat.raw.pbulk <- SumAcrossClusters(count.mat, cnames.keep.lst)
  dat.raw.pbulk <- do.call(cbind, dat.raw.pbulk)
  # normalize
  dat.raw.pbulk <- sweep(dat.raw.pbulk, MARGIN = 2, STATS = colSums(dat.raw.pbulk), FUN = "/", check.margin = TRUE) * jscale + 1
})

# save to output?
# saveRDS(dat.raw.pbulk.lst, file = paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pseudobulk_high_bins_matrix/dat_raw_four_marks_high_bins_pseudobulk.", Sys.Date(), ".rds"))



# Calculate variance for each gene attribute to cluster -------------------

# dat.raw.pbulk.lst <- dat.raw.pbulk.lst.lst$DE_bins_all_marks_top_6085_dists_to_TSS.annot_table

head(dat.raw.pbulk.lst[[1]])
head(log2(dat.raw.pbulk.lst[[1]]))

jmat <- log2(dat.raw.pbulk.lst[[1]])
jcheck <- sweep(jmat, MARGIN = 1, STATS = rowMeans(jmat), FUN = "-") 

dat.vars.lst <- lapply(dat.raw.pbulk.lst, function(jdat){
  jdat.log2 <- log2(jdat)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
  jcheck <- jcheck ^ 2
})


dat.deviations.lst <- lapply(dat.raw.pbulk.lst, function(jdat){
  jdat.log2 <- log2(jdat)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
})


dat.deviations.relativetohspcs.lst <- lapply(dat.deviations.lst, function(jdat){
  jdat.rel <- sweep(jdat, MARGIN = 1, STATS = jdat[, "HSPCs"], FUN = "-")
  jcols.keep <- colnames(jdat.rel) != "HSPCs"
  return(jdat.rel[, jcols.keep])
})


min.log2fc <- 0
dat.deviations.mean.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.deviations.relativetohspcs.lst[[jmark]]
  rnames.keep <- apply(abs(jdat), MARGIN = 1, FUN = function(jrow) max(jrow) > min.log2fc)
  jdat <- jdat[rnames.keep, ]
  jdat.long <- data.frame(bin = rownames(jdat), mean.dev = rowMeans(jdat))
  jdat.long$mark <- jmark
  return(jdat.long)
})


dat.deviations.mean.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.deviations.relativetohspcs.lst[[jmark]]
  rnames.keep <- apply(abs(jdat), MARGIN = 1, FUN = function(jrow) max(jrow) > min.log2fc)
  jdat <- jdat[rnames.keep, ]
  jdat.long <- data.frame(bin = rownames(jdat), mean.dev = rowMeans(jdat))
  jdat.long$mark <- jmark
  return(jdat.long)
})

# deviations from 



# Plot mean deviations from HSPCs -----------------------------------------

ggplot(dat.deviations.mean.lst %>% bind_rows(), aes(x = mean.dev, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  facet_wrap(~mark, ncol = 1) + 
  geom_vline(xintercept = c(-2.5, 0, 2.5), linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load dynamic bins  ------------------------------------------------------

infs.de.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.de.bins.lst <- lapply(infs.de.lst, function(jinf){
  fread(jinf)
})

bins.de.lst <- lapply(dat.de.bins.lst, function(jdat) jdat$CoordOriginal)

dat.deviations.mean.filt.long <- lapply(jmarks, function(jmark){
  jdat <- dat.deviations.mean.lst[[jmark]]
  bins.keep <- bins.de.lst[[jmark]]
  jdat.filt <- subset(jdat, bin %in% bins.keep)
})  %>%
  bind_rows()
dat.deviations.mean.filt.long$mark <- factor(dat.deviations.mean.filt.long$mark, levels = jmarks)


ggplot(dat.deviations.mean.filt.long, aes(x = mean.dev, fill = mark)) + 
  geom_density(alpha = 0.25, fill = "grey25") + 
  xlab("Deviation from HSPCs averaged across cell types [log2FC]") + 
  theme_bw() + 
  facet_wrap(~mark, nrow = 1) + 
  geom_vline(xintercept = c(0), linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")





# Calculate variance of each gene  ----------------------------------------

dat.vars.long.lst <- lapply(jmarks, function(jmark){
  dat.vars <- dat.vars.lst[[jmark]]
  melt(dat.vars) %>%
    dplyr::rename(bin = Var1, 
                  celltype = Var2) %>%
    group_by(celltype) %>%
    summarise(vartotal = sum(value)) %>%
    mutate(mark = jmark) %>%
    group_by(mark) %>%
    mutate(varfrac = vartotal / sum(vartotal)) %>%
    rowwise() %>%
    mutate(lineage = ctypes[[celltype]],
           clustercol = cluster2col[[as.character(celltype)]])
})


# Do bootstraps -----------------------------------------------------------

set.seed(0)
nbootstraps.vec <- seq(1000); names(nbootstraps.vec) <- nbootstraps.vec
nsamps <- 5000

# jmarksub <- jmarks[1]
dat.vars.long.subsamp.lst.lst <- lapply(nbootstraps.vec, function(i){
  if (i %% 100 == 0){
    print(i)
  }
  dat.vars.long.subsamp.lst <- lapply(jmarks, function(jmark){
    # print(jmark)
    dat.vars <- dat.vars.lst[[jmark]]
    # subsamp
    rows.subsamp <- sample(x = seq(nrow(dat.vars)), size = nsamps, replace = FALSE)
    dat.vars <- dat.vars[rows.subsamp, ]
    
    dat.vars.sum.long <- melt(dat.vars) %>%
      dplyr::rename(bin = Var1, 
                    celltype = Var2) %>%
      group_by(celltype) %>%
      summarise(vartotal = sum(value)) %>%
      mutate(mark = jmark) %>%
      group_by(mark) %>%
      mutate(varfrac = vartotal / sum(vartotal)) %>%
      rowwise() %>%
      mutate(varnorm = vartotal / nsamps, 
             lineage = ctypes[[celltype]],
             clustercol = cluster2col[[as.character(celltype)]])
    dat.vars.sum.long$i <- i
    return(dat.vars.sum.long)
  })
})

# plot bootstraps

prob.low <- 0.01
prob.high <- 0.99
dat.vars.long.subsamp.long.sum.lst <- lapply(jmarks, function(jmark){
  dat.vars.long.subsamp.long.sum <- lapply(dat.vars.long.subsamp.lst.lst, function(x){
    x[[jmark]]
  }) %>%
    bind_rows() %>%
    group_by(celltype, mark) %>%
    summarise(jmedian = median(varnorm),
              jlower = quantile(varnorm, probs = prob.low),
              jupper = quantile(varnorm, probs = prob.high))
})






# Normalize and plot ------------------------------------------------------

nbins.lst <- lapply(jmarks, function(jmark){
  nrow(dat.raw.pbulk.lst[[jmark]])
})

# add error bars 
dat.vars.long.annot.lst <- lapply(jmarks, function(jmark){
  jmerge <- left_join(dat.vars.long.lst[[jmark]], dat.vars.long.subsamp.long.sum.lst[[jmark]], by = c("celltype", "mark"))
  jmerge$nbins <- nbins.lst[[jmark]]
  jmerge$varnorm <- jmerge$vartotal / jmerge$nbins
  return(jmerge)
})

dat.vars.long <- dat.vars.long.annot.lst %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(devnorm = sqrt(varnorm),
         fcnorm = 2^devnorm)

# # normalize var
# dat.vars.long <- dat.vars.long %>%
#   group_by(mark) %>%
#   mutate(varfrac = vartotal / sum(vartotal))

head(dat.vars.long)

subset(dat.vars.long, mark == "H3K4me1")

ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = varnorm, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  ggtitle(paste("Error bars:", prob.low, prob.high)) + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  geom_errorbar(mapping = aes(ymin = jlower, ymax = jupper), color = "grey25", width = 0.5) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = devnorm, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  ggtitle(paste("Error bars:", prob.low, prob.high)) + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = fcnorm, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  ggtitle(paste("Error bars:", prob.low, prob.high)) + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))




ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.vars.long, aes(x = forcats::fct_reorder(.f = celltype, .x = varfrac, .fun = median, .desc = TRUE), y = varfrac, fill = clustercol)) + 
  geom_col() + 
  scale_fill_identity() + 
  facet_wrap(~mark, scales = "free_x") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


for (jmark in jmarks){
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
    geom_col() + 
    ylim(c(0, jymax)) + 
    scale_fill_identity() + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}


# ylim(c(0, jymaxnorm)) + 
for (jmark in jmarks){
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = varnorm, fill = clustercol)) + 
    geom_col() + 
    geom_errorbar(mapping = aes(ymin = jlower, ymax = jupper), color = "grey25", width = 0.5) + 
    scale_fill_identity() + 
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, jymaxnorm)) + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.vars.long %>% filter(mark == jmark), aes(x = forcats::fct_reorder(.f = celltype, .x = varfrac, .fun = median, .desc = TRUE), y = varfrac, fill = clustercol)) + 
    geom_col() + 
    ylim(c(0, jymaxfrac)) + 
    scale_fill_identity() + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

# split variance by pos or neg


dat.diffs.lst <- lapply(dat.raw.pbulk.lst, function(jdat){
  # jdat <- jout$dat.impute.pbulk
  # jdat <- jout$dat.raw.pbulk
  jdat.log2 <- log2(jdat)
  jcheck <- sweep(jdat.log2, MARGIN = 1, STATS = rowMeans(jdat.log2), FUN = "-") 
  # jcheck <- jcheck ^ 2
})

# calculate variance, break down into pos or neg
dat.vars.posneg.long.lst <- lapply(jmarks, function(jmark){
  dat.diffs <- dat.diffs.lst[[jmark]]
  dat.diffs.long <- melt(dat.diffs)
  colnames(dat.diffs.long) <- c("rname", "celltype", "log2diff")
  dat.diffs.sum.long <- dat.diffs.long %>%
    rowwise() %>%
    mutate(is.pos = log2diff > 0) %>%
    group_by(celltype, is.pos) %>%
    summarise(vartotal = sum(log2diff ^ 2),
              nbins = length(log2diff)) %>%
    rowwise() %>%
    mutate(varnorm = vartotal / nbins) %>%
    mutate(lineage = ctypes[[celltype]],
           clustercol = cluster2col[[as.character(celltype)]])
  dat.diffs.sum.long$mark <- jmark
  return(dat.diffs.sum.long)
})



for (jmark in jmarks){
  m <- ggplot(dat.vars.long.lst[[jmark]], aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = clustercol)) + 
    geom_col() + 
    ylim(c(0, jymax)) + 
    scale_fill_identity() + 
    xlab("") + 
    ylab("Total Variance") + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.vars.posneg.long.lst[[jmark]], aes(x = forcats::fct_reorder(.f = celltype, .x = vartotal, .fun = median, .desc = TRUE), y = vartotal, fill = is.pos)) + 
    # geom_col(position = "stack") + 
    geom_col() + 
    ylim(c(0, jymax)) + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    xlab("") + 
    ylab("Total Variance") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")
  print(m)
}


for (jmark in jmarks){
  m <- ggplot(dat.vars.posneg.long.lst[[jmark]], aes(x = forcats::fct_reorder(.f = celltype, .x = varnorm, .fun = median, .desc = TRUE), y = varnorm, fill = is.pos)) + 
    # geom_col(position = "stack") + 
    geom_col(position = "dodge") + 
    ylim(c(0, jymaxnorm * 1.5)) + 
    facet_wrap(~mark, scales = "free_x") + 
    theme_bw() + 
    xlab("") + 
    ylab("Variance per bin, split by pos or not") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")
  print(m)
}



# make heatmaps 
# 
# # jmark <- "H3K9me3"
# jmethods <- c("complete", "ward.D2")
# # jmethod <- "ward.D2"
# # jmethod <- "complete"
# for (jmethod in jmethods){
#   
#   for (jmark in jmarks){
#     print(paste("Making heatmap for", jmark))
#     jmat2 <- log2(dat.raw.pbulk.lst[[jmark]])
#     print(dim(jmat2))
#     Ngenes <- nrow(jmat2)
#     cnames.color <- sapply(colnames(jmat2), function(x) AssignHash(x, cname2color, null.fill = x))
#     heatmap3::heatmap3(jmat2, 
#                        Rowv = TRUE, Colv = TRUE, scale = "row",  revC = TRUE, 
#                        main = paste("peaks", jmark, jmethod, Ngenes), margins = c(5, 8), cexRow = 0.5, method = jmethod, ColSideColors = cnames.color, col = colpalette)
#   }
# }

if (make.plots){
  dev.off()
}
