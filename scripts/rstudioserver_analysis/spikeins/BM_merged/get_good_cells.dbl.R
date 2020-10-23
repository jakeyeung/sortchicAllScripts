# Jake Yeung
# Date of Creation: 2020-10-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/get_good_cells.R
# Filter good cells

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"
# indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/countTablesAndRZr1only_ByChromo.NewFilters")
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams/dbl_stains/countTablesAndRZr1only_ByChromo.NewFilters.blfix")
indir.others <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams/dbl_stains/countTablesAndRZr1only_ByChromo.NewFilters")
assertthat::assert_that(dir.exists(indir))


jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")
jmarks <- c("H3K9me3-H3K4me1"); names(jmarks) <- jmarks

# Load everything ---------------------------------------------------------

infs.bins <- list.files(path = indir, pattern = "*binsize_50000.csv", full.names = TRUE)
names(infs.bins) <- sapply(infs.bins, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

infs.rz <- list.files(path = indir.others, pattern = "*countTable.RZ.csv", full.names = TRUE)
names(infs.rz) <- sapply(infs.rz, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

infs.chromo <- list.files(path = indir.others, pattern = "*ByChromo.WithSpikeIns.NoChromo.csv", full.names = TRUE)
# names(infs.chromo) <- gsub(pattern = "H4K4me1", replacement = "H3K4me1", sapply(infs.chromo, basename))
names(infs.chromo) <- sapply(infs.chromo, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

assertthat::assert_that(identical(names(infs.rz), names(infs.chromo)))

names.all <- names(infs.rz)

dats.rzs <- lapply(infs.rz, function(inf.rz){
  print(inf.rz)
  dat.rz <- ReadLH.SummarizeTA(inf.rz)
}) %>%
  bind_rows()

dats.chromos <- lapply(infs.chromo, function(inf.chromo){
  print(inf.chromo)
  dat.chromos <- GetChromoCounts(inf.chromo, spikeinchromo = jspikeinchromo, chromos.keep = jchromos) %>%
    filter(chromo == "1")
}) %>%
  bind_rows()

dats.rzchromo <- left_join(dats.rzs, dats.chromos, by = "samp")




dats.rzchromo <- dats.rzchromo %>%
  rowwise() %>%
  mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
         colcoord = AddPlateCoordinates(samp)$colcoord,
         is.empty = rowcoord <= 8 & colcoord == 1, 
         mark = jmarks[[1]])


# Get good cells  ---------------------------------------------------------


# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl.blfix"
dir.create(outdir)
outpdf <- file.path(outdir, "qc_plots_all_marks_with_spikeins.blfix.pdf")
# assertthat::assert_that(!file.exists(outpdf))




# Get good cells  ---------------------------------------------------------

# remove 15% of cells left and 5% right from each plate? 

dats.rzchromo.annot <- dats.rzchromo %>%
  group_by(mark) %>%
  mutate(chromocounts.rnk = rank(chromocounts) / length(chromocounts), 
         l2r = log2(chromocounts / spikeincounts), 
         l2r.rnk = rank(l2r) / length(l2r)) %>%
  rowwise() %>%
  # mutate(is.good = !is.empty & chromocounts.rnk >= 0.15 & chromocounts.rnk <= 0.05 & l2r.rnk >= 0.15 & l2r.rnk <= 0.05)
  mutate(is.good = !is.empty & TA.frac > 0.5 & chromocounts.rnk >= 0.15 & chromocounts.rnk <= 0.95 & l2r.rnk >= 0.15 & l2r.rnk <= 0.95)




pdf(outpdf, useDingbats = FALSE)
mlst <- lapply(jmarks, function(jmark){
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "all")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark & is.good), aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "filt")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), fill = is.empty)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "all")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark & is.good), aes(x = log2(chromocounts / spikeincounts), fill = is.empty)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "filt")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark), aes(x = log10(chromocounts), y = TA.frac, color = is.empty)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "all")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark & is.good), aes(x = log10(chromocounts), y = TA.frac, color = is.empty)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "filt")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark & is.good), aes(x = log10(chromocounts), fill = is.empty)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "all")
  print(m)
  
  m <- ggplot(dats.rzchromo.annot %>% filter(mark == jmark & is.good), aes(x = log10(chromocounts), fill = is.empty)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark, "filt")
  print(m)
  
})
dev.off()

# Load mats ---------------------------------------------------------------


mats <- lapply(infs.bins, function(inf){
  print(inf)
  ReadMatSlideWinFormat(inf)
})

names.all <- names(mats)

rnames.all <- Reduce(f = intersect, x = lapply(names.all, function(x) rownames(mats[[x]])))

# Merge by mark  ----------------------------------------------------------

names.bymark <- list(names.all)
names(names.bymark) <- jmarks

rnames.bymark <- list(rnames.all)
names(rnames.bymark) <- jmarks

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

# keep mouse chromos only
rnames.bymark.filt <- lapply(rnames.bymark, function(rnames){
  chromos.vec <- sapply(rnames, GetChromo)
  # keep only that are in jchromos
  keepi <- chromos.vec %in% jchromos
  rnames.keep <- rnames[keepi]
  return(rnames.keep)
})

cells.bymark <- lapply(jmarks, function(jmark){
  subset(dats.rzchromo.annot, mark == jmark & is.good)$samp
})

mats.bymark <- lapply(jmarks, function(jmark){
  print(jmark)
  jnames <- names.bymark[[jmark]]
  rnames.keep <- rnames.bymark.filt[[jmark]]
  cells.keep <- cells.bymark[[jmark]]
  mat.merged <- lapply(jnames, function(jname){
    print(jname)
    mat <- mats[[jname]]
    return(mat[rnames.keep, ])
  })
  mat.merged <- do.call(cbind, mat.merged)
  print(dim(mat.merged))
  # print(mat.merged)
  cols.keep <- colnames(mat.merged) %in% cells.keep
  return(mat.merged[, cols.keep])
})


# Write output ------------------------------------------------------------


for(jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat.", jmark, ".filt_0.15_0.95_counts_and_l2r.blfix.rds"))
  # assertthat::assert_that(!file.exists(outf))
  
  jmat <- mats.bymark[[jmark]]
  print(dim(jmat))
  saveRDS(object = jmat, file = outf)
}


# Write spikein -----------------------------------------------------------

outf.spikein <- file.path(outdir, "spikein_info_BM_round2_all.blfix.txt")
assertthat::assert_that(!file.exists(outf.spikein))
fwrite(dats.rzchromo, file = outf.spikein)







