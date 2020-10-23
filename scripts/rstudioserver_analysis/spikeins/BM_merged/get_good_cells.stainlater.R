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
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams/BM_stain_later/countTablesAndRZr1only_ByChromo.NewFilters.blfix")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.stainlater"
dir.create(outdir)
outpdf <- file.path(outdir, "qc_plots_all_marks_with_spikeins.pdf")

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

# Load everything ---------------------------------------------------------

infs.bins <- list.files(path = indir, pattern = "*binsize_50000.csv", full.names = TRUE)
names(infs.bins) <- sapply(infs.bins, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

infs.rz <- list.files(path = indir, pattern = "*countTable.RZ.csv", full.names = TRUE)
names(infs.rz) <- sapply(infs.rz, function(x) strsplit(basename(x), split = "\\.")[[1]][[1]])

infs.chromo <- list.files(path = indir, pattern = "*ByChromo.WithSpikeIns.NoChromo.csv", full.names = TRUE)
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



GetMarkFromSamp2 <- function(samp){
  if (grepl("K4me1", samp)){
    return("H3K4me1")
  } else if (grepl("K4me3", samp)){
    return("H3K4me3")
  } else if (grepl("K27me3", samp)){
    return("H3K27me3")
  } else if (grepl("K9me3", samp)){
    return("H3K9me3")
  }
}

dats.rzchromo <- dats.rzchromo %>%
  rowwise() %>%
  mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord,
         colcoord = AddPlateCoordinates(samp)$colcoord,
         is.empty = rowcoord <= 8 & colcoord == 1, 
         mark = GetMarkFromSamp2(samp))


# Get good cells  ---------------------------------------------------------

jmarks <- c("H3K27me3"); names(jmarks) <- jmarks


pdf(outpdf, useDingbats = FALSE)
mlst <- lapply(jmarks, function(jmark){
  m <- ggplot(dats.rzchromo %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), y = TA.frac, color = is.empty)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark)
  print(m)
    
  m <- ggplot(dats.rzchromo %>% filter(mark == jmark), aes(x = log2(chromocounts / spikeincounts), fill = is.empty)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(dats.rzchromo %>% filter(mark == jmark), aes(x = log10(chromocounts), y = TA.frac, color = is.empty)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark)
  print(m)
  
  m <- ggplot(dats.rzchromo %>% filter(mark == jmark), aes(x = log10(chromocounts), fill = is.empty)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    coord_cartesian(ylim = c(0, 1)) + 
    facet_wrap(~experi) + 
    ggtitle(jmark)
  print(m)
})
dev.off()


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



# Load mats ---------------------------------------------------------------


mats <- lapply(infs.bins, function(inf){
  print(inf)
  ReadMatSlideWinFormat(inf)
})

names.k27me3 <- grep(pattern = "K27me3", names(mats), value = TRUE)
names(names.k27me3) <- names.k27me3

rnames.k27me3 <- Reduce(f = intersect, x = lapply(names.k27me3, function(x) rownames(mats[[x]])))



# Merge by mark  ----------------------------------------------------------

names.bymark <- list(names.k27me3)
names(names.bymark) <- jmarks

rnames.bymark <- list(rnames.k27me3)
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
  outf <- file.path(outdir, paste0("count_mat.", jmark, ".filt_0.15_0.95_counts_and_l2r.rds"))
  jmat <- mats.bymark[[jmark]]
  print(dim(jmat))
  saveRDS(object = jmat, file = outf)
}


# Write spikein -----------------------------------------------------------

outf.spikein <- file.path(outdir, "spikein_info_BM_round2_all.txt")
fwrite(dats.rzchromo, file = outf.spikein)







