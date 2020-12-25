# Jake Yeung
# Date of Creation: 2020-10-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/add_hoescht_to_cellcycle_K562.R
# For AvO



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

qvalcutoff <- 3


outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/glmpca_spikeins_hoescht.K562_G1_G2/glmpca_spikeins_hoescht_summaries.K562_G1_G2.", Sys.Date(), ".PeakCounts.qval_", qvalcutoff, ".pdf")
outrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/glmpca_spikeins_hoescht.K562_G1_G2/glmpca_spikeins_hoescht_summaries.K562_G1_G2.", Sys.Date(), ".PeakCounts.qval_", qvalcutoff, ".RData")

assertthat::assert_that(!file.exists(outrdata))
assertthat::assert_that(!file.exists(outpdf))


# Get Hoesch information --------------------------------------------------


# Load spikeins -----------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.rds"
outf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.WithCellCycleLabel.rds"
dat.spikeins <- readRDS(inf.spikeins)
dat.spikeins <- AddCellCycleLabel.bydat(dat.spikeins)

saveRDS(dat.spikeins, file = outf.spikeins)

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


# Save to output ----------------------------------------------------------

outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/K562_RData_objs/K562_cellcycle_with_hoescht_4_marks.rds"
saveRDS(dat.glmpca.lst, file = outf)
