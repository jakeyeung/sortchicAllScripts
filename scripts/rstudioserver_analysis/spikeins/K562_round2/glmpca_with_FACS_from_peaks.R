# Jake Yeung
# Date of Creation: 2020-10-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/glmpca_with_FACS.R
# 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(irlba)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jmarks.hoescht.prefix <- c("4me1", "4me3", "27me3", "9me3")
hoescht.suffix <- c(".*.4.*index$")
jmarks.hoescht <- paste(jmarks.hoescht.prefix, hoescht.suffix, sep = "")
names(jmarks.hoescht) <- jmarks

make.plots <- TRUE


# Get Hoesch information --------------------------------------------------


# Load spikeins -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.rds"
dat.spikeins <- readRDS(inf.spikeins)

dat.spikeins <- AddCellCycleLabel.bydat(dat.spikeins)

# Load LDA  ---------------------------------------------------------------

# jmark <- jmarks[[3]]
jsuffix <- "cellcyclefilt"
# jsuffix <- "AllMerged"
# jsuffix <- "G1filt"
jtopn <- 5000
jmethod <- "avagrad"
jminibatch <- "stochastic"
# jtol <- "1e-8"
jtol <- "1e-4"
jiter <- 500

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

if (jsuffix == "AllMerged"){
  allmerged <- TRUE
} else {
  allmerged <- FALSE
}

# if (!allmerged){
#   outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562.round2.WithFACS.peaks"
# } else {
#   outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562.round2.AllMerged.WithFACS.peaks"
# }
# dir.create(outdir)


inf.hoescht <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.hoesch_staining/hoesch_on_K562_plates.rds"


dat.glmpca.lst <- lapply(jmarks, function(jmark){
  # Load glmpca -------------------------------------------------------------
  
  # inf.glmpca <- paste0("/home/jyeu  ng/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2/K562_count_tables_50000.", jmark, ".", jsuffix, ".glmpcaout.penalty_1.maxiter_10000.RData")
  if (!allmerged){
    inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2.genes_filt/K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".glmpcaout.penalty_1.maxiter_", jiter, ".", jminibatch, ".", jmethod, ".tol_", jtol, ".RData")
  } else {
    inf.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2.merged.genes_filt/K562_count_tables_50000.", jmark, ".", jsuffix, ".topn_", jtopn, ".glmpcaout.penalty_1.maxiter_", jiter, ".", jminibatch, ".", jmethod, ".tol_", jtol, ".RData")
  }
  assertthat::assert_that(file.exists(inf.glmpca))
  load(inf.glmpca, v=T)
  
  # Plot GLMPCA -------------------------------------------------------------
  
  dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings)
  dat.umap.glmpca$dim1 <- glmpcaout$factors$dim1
  dat.umap.glmpca$dim2 <- glmpcaout$factors$dim2
  dat.umap.glmpca$mark <- jmark
  
  # add FACS
  
  
  dat.umap.glmpca <- dat.umap.glmpca %>%
    rowwise() %>%
    mutate(experi2 = strsplit(cell, "_")[[1]][[1]])  %>%
    left_join(., dat.spikeins, by = c("cell" = "samp"))  # add meta data
  
  
  
  experi.str <- unique(dat.umap.glmpca$experi2)
  assertthat::assert_that(length(experi.str) == 1)
  
  jmark.hoescht <- jmarks.hoescht[[jmark]]
  print(paste(jmark, ":", jmark.hoescht))
  
  dat.hoescht.all <- readRDS(inf.hoescht)
  dat.hoescht <- subset(dat.hoescht.all, grepl(jmark.hoescht, x = experi))
  
  dat.hoescht <- dat.hoescht %>%
    rowwise() %>%
    mutate(experi2 = experi.str)
  
  # rename experi to match glmpca
  dat.hoescht <- dat.hoescht %>%
    rowwise() %>%
    mutate(experi2 = experi.str) %>%
    ungroup() %>%
    mutate(hoesch.log = hoesch,
           hoesch.scale = (hoesch.log - min(hoesch.log)) / (max(hoesch.log) - min(hoesch.log)))
  
  dat <- left_join(dat.umap.glmpca, dat.hoescht, by = c("rowcoord" = "row.indx", "colcoord" = "col.indx", "experi2" = "experi2"))
  return(dat)
})


# dat.glmpca.merged <- left_join(dat.glmpca.lst %>% bind_rows(), dat.spikeins, by = c("cell" = "samp"))
ggplot(dat.glmpca.lst[["H3K4me1"]], aes(x = umap1, y = umap2, color = cellcycle.str)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.glmpca.lst[["H3K4me1"]], aes(x = umap1, y = umap2, color = hoesch)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jmark <- jmarks[[1]]
m <- ggplot(dat.glmpca.lst[[jmark]], aes(x = hoesch.scale, y = log2(chromocounts / spikeincounts), color = cellcycle.str)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle(jmark) + 
  facet_wrap(~experi2) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)


jslopes <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.lst[[jmark]]
  jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ hoesch.scale, data = jsub)
  jslope <- summary(jfit)$coefficients["hoesch.scale", "Estimate"]
  jslope.dat <- data.frame(slope = jslope, mark = jmark, stringsAsFactors = FALSE)
  return(jslope.dat)
}) 

# 
# jslopes.discrete <- lapply(jmarks, function(jmark){
#   jsub <- dat.glmpca.lst[[jmark]]
#   jfit <- lm(formula = log2(peakcounts / spikeincounts) ~ cellcycle.str, data = jsub)
#   jslope <- summary(jfit)$coefficients["cellcycle.str2_G2/M", "Estimate"]
#   jslope.dat <- data.frame(slope = jslope, mark = jmark, stringsAsFactors = FALSE)
#   return(jslope.dat)
# }) 
# 
# 

if (make.plots){
  outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/glmpca_spikeins_hoescht.K562_G1_G2/glmpca_spikeins_hoescht_summaries.K562_G1_G2.", Sys.Date(), ".PeakCounts.pdf")
  pdf(outpdf, file = outpdf)
}

mlst <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.lst[[jmark]]
  jsub$ltr <- log2(jsub$chromocounts / jsub$spikeincounts)
  jsub$ltr.wins <- DescTools::Winsorize(jsub$ltr, probs = c(0.01, 0.99))
  
  m.umap <- ggplot(jsub, aes(x = umap1, y = umap2, color = ltr.wins)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "log2(cuts/spikeins), winsorized") + 
    facet_wrap(~experi2) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.umap)
  
  m.dim <- ggplot(jsub, aes(x = dim1, y = dim2, color = ltr.wins)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "log2(cuts/spikeins), winsorized") + 
    facet_wrap(~experi2) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.dim)
 
  
  jslope <- jslopes[[jmark]]$slope
  m <- ggplot(jsub, aes(x = hoesch.scale, y = log2(chromocounts / spikeincounts), color = cellcycle.str)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~experi2) + 
    ggtitle(paste("Slope:", signif(jslope, digits = 2))) + 
    geom_smooth(mapping = aes(x = hoesch.scale, y = log2(chromocounts / spikeincounts)), method = "lm", se = FALSE, inherit.aes = FALSE, color = 'black', alpha = 0.75) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(jsub, aes(x = hoesch.scale, y = ltr.wins, color = cellcycle.str)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~experi2) + 
    ylab("log2(cuts/spikeins), winsorized") + 
    geom_smooth(mapping = aes(x = hoesch.scale, y = log2(chromocounts / spikeincounts)), method = "lm", se = FALSE, inherit.aes = FALSE, color = 'black', alpha = 0.75) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m.small <- ggplot(jsub, aes(x = hoesch.scale, y = log2(chromocounts / spikeincounts), color = cellcycle.str)) + 
    geom_point(alpha = 0.25) + 
    theme_bw(8) + 
    ggtitle(jmark) + 
    facet_wrap(~experi2) + 
    ggtitle(paste("Slope:", signif(jslope, digits = 2))) + 
    geom_smooth(mapping = aes(x = hoesch.scale, y = log2(chromocounts / spikeincounts)), method = "lm", se = FALSE, inherit.aes = FALSE, color = 'black', alpha = 0.75) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m.small)
})

JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)





# Calculate new chromo counts  --------------------------------------------

good.cells.lst <- lapply(dat.glmpca.lst, function(jdat){
  jdat$cell
})

hubprefix <- "/home/jyeung/hub_oudenaarden"
indir.counts <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBed.NewFilters")

mats.counts <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir.counts, paste0("K562_CellCycleSorted_", jmark, ".merged.sorted.tagged.bam.countTable.bedfile.csv"))
  ReadMatTSSFormat(inf.tmp)
})

# get good cells 
mats.counts.filt <- lapply(jmarks, function(jmark){
  jmat <- mats.counts[[jmark]]
  cols.keep <- colnames(jmat) %in% good.cells.lst[[jmark]]
  assertthat::assert_that(length(which(cols.keep)) > 0)
  jmat[, cols.keep]
})

chromo.counts.peaks <- lapply(mats.counts.filt, function(jmat){
  data.frame(cell = colnames(jmat), peakcounts = colSums(jmat), stringsAsFactors = FALSE)
})

# add to dat.glmpca

dat.glmpca.peak.lst <- lapply(jmarks, function(jmark){
  left_join(dat.glmpca.lst[[jmark]], chromo.counts.peaks[[jmark]])
}) 

m.lst <- lapply(dat.glmpca.peak.lst, function(jdat){
  m <- ggplot(jdat, aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

JFuncs::multiplot(m.lst$H3K4me1, m.lst$H3K4me3, m.lst$H3K27me3, m.lst$H3K9me3, cols = 4)


# Refit -------------------------------------------------------------------

dat.glmpca.peak.lst <- lapply(dat.glmpca.peak.lst, function(jsub){
  jsub$ltr <- log2(jsub$peakcounts / jsub$spikeincounts)
  jsub$ltr.wins <- DescTools::Winsorize(jsub$ltr, probs = c(0.01, 0.99))
  return(jsub)
})

jslopes.peak <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.peak.lst[[jmark]]
  jfit <- lm(formula = log2(peakcounts / spikeincounts) ~ hoesch.scale, data = jsub)
  jslope <- summary(jfit)$coefficients["hoesch.scale", "Estimate"]
  jslope.dat <- data.frame(slope = jslope, mark = jmark, stringsAsFactors = FALSE)
  return(jslope.dat)
}) 

jslopes.peak.winsorized <- lapply(jmarks, function(jmark){
  jsub <- dat.glmpca.peak.lst[[jmark]]
  jfit <- lm(formula = ltr.wins ~ hoesch.scale, data = jsub)
  jslope <- summary(jfit)$coefficients["hoesch.scale", "Estimate"]
  jslope.dat <- data.frame(slope = jslope, mark = jmark, stringsAsFactors = FALSE)
  return(jslope.dat)
}) 


mlst.peak <- lapply(jmarks, function(jmark){
  print(jmark)
  jsub <- dat.glmpca.peak.lst[[jmark]]
  # remove two outliers in K27me3?
  if (jmark == "H3K27me3"){
    jsub <- subset(jsub, log2(peakcounts/spikeincounts) > 2)
  }
  jsub$ltr <- log2(jsub$peakcounts / jsub$spikeincounts)
  jsub$ltr.wins <- DescTools::Winsorize(jsub$ltr, probs = c(0.01, 0.99))
  
  m.umap <- ggplot(jsub, aes(x = umap1, y = umap2, color = ltr.wins)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "log2(cuts/spikeins), winsorized") + 
    facet_wrap(~experi2) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.umap)
  
  m.dim <- ggplot(jsub, aes(x = dim1, y = dim2, color = ltr.wins)) + 
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "log2(cuts/spikeins), winsorized") + 
    facet_wrap(~experi2) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.dim)
  
  
  jslope <- jslopes.peak[[jmark]]$slope
  m <- ggplot(jsub, aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts), color = cellcycle.str)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~experi2) + 
    ggtitle(paste("Slope:", signif(jslope, digits = 2))) + 
    geom_smooth(mapping = aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts)), method = "lm", se = FALSE, inherit.aes = FALSE, color = 'black', alpha = 0.75) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(jsub, aes(x = hoesch.scale, y = ltr.wins, color = cellcycle.str)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~experi2) + 
    ylab("log2(cuts/spikeins), winsorized") + 
    geom_smooth(mapping = aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts)), method = "lm", se = FALSE, inherit.aes = FALSE, color = 'black', alpha = 0.75) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m.small <- ggplot(jsub, aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts), color = cellcycle.str)) + 
    geom_point(alpha = 0.25) + 
    theme_bw(8) + 
    ggtitle(jmark) + 
    facet_wrap(~experi2) + 
    ggtitle(paste("Slope:", signif(jslope, digits = 2))) + 
    geom_smooth(mapping = aes(x = hoesch.scale, y = log2(peakcounts / spikeincounts)), method = "lm", se = FALSE, inherit.aes = FALSE, color = 'black', alpha = 0.75) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m.small)
})

JFuncs::multiplot(mlst.peak[[1]], mlst.peak[[2]], mlst.peak[[3]], mlst.peak[[4]], cols = 4)



if (make.plots){
  dev.off()
}


# Save object for downstream ----------------------------------------------


outrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/glmpca_spikeins_hoescht.K562_G1_G2/glmpca_spikeins_hoescht_summaries.K562_G1_G2.", Sys.Date(), ".PeakCounts.RData")
save(jslopes.peak, dat.glmpca.peak.lst, mats.counts.filt, file = outrdata)


