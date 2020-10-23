# Jake Yeung
# Date of Creation: 2020-10-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/H3K27me3_merged/2-LDA_downstream.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)




AnnotateSortFromLayout <- function(plate, rowcoord, colcoord){
  assertthat::assert_that(is.numeric(plate))
  assertthat::assert_that(is.numeric(rowcoord))
  assertthat::assert_that(is.numeric(colcoord))
  if (rowcoord >= 1 & rowcoord <= 8 & colcoord == 1){
    print("Empty cell")
    warning("Empty cell")
    ctype <- "Empty"
    return(ctype)
  }
  if (plate >= 1 & plate <= 7){
    if (colcoord >= 1 & colcoord <= 11){
      ctype <- "Unenriched"
    } else if (colcoord >= 12 & colcoord <= 18){
      ctype <- "LinNeg"
    } else if (colcoord >= 19 & colcoord <= 24){
      ctype <- "LSK"
    }
    
  } else if (plate >= 8 & plate <= 13){
    if (colcoord >= 1 & colcoord <= 12){
      ctype <- "Unenriched"
    } else if (colcoord >= 13 & colcoord <= 24){
      ctype <- "LinNeg"
    }
  } else {
    print(paste("Unknown plate:", plate))
    ctype <- "Unknown"
  }
  return(ctype)
}



jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

var.cutoff <- 1.5
fname <- "count_mat_H3K27me3_l2r_filt.2020-10-03.minl2r_0"


# Load data  --------------------------------------------------------------

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.merged_across_runs/varfilt"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.round2/VAN5046_VAN5230_BM"
dir.create(outdir)
# outf <- file.path(outdir, paste0(fname, ".varfilt_", var.cutoff, ".rds"))
outpdf <- file.path(outdir, paste0(fname, ".varfilt_", var.cutoff, ".downstream.pdf"))

pdf(outpdf, useDingbats = FALSE)

# pdf(file = outpdf, useDingbats = FALSE)
# 
inf.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046_VAN5230_BM/lda_outputs.", fname, ".K-30.binarize.FALSE/ldaOut.", fname, ".K-30.Robj")
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result, jsep = "")

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# get variance
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.merge <- left_join(dat.umap, dat.var)

ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Add whether it is unenriched, linneg or not -----------------------------

dat.merge.annot <- dat.merge %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))

ggplot(dat.merge.annot, aes(x = umap1, y = umap2, color = stype)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~experi) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 

# Losad spikeins  ---------------------------------------------------------

inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.merged_across_runs/spikeins_dat_H3K27me3_merged.txt"
dat.spikeins <- fread(inf.spikeins)

dat.merge.annot2 <- left_join(dat.merge.annot, dat.spikeins)
dat.merge.annot2$stype <- factor(dat.merge.annot2$stype, levels = c("LSK", "LinNeg", "Unenriched"))

ggplot(dat.merge.annot2, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

ggplot(dat.merge.annot2, aes(x = log2(chromocounts / spikeincounts), y = cell.var.within.sum.norm)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

# Check spikeins by plate ---------------------------------------------------------


ggplot(dat.merge.annot2, aes(x = stype, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_jitter(mapping = aes(color = as.character(plate)), width = 0.1) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


ggplot(dat.merge.annot2, aes(x = stype, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("") + 
  facet_wrap(~plate)

ggplot(dat.merge.annot2, aes(x = stype, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("") + 
  facet_wrap(~plate, nrow = 1)
  
ggplot(dat.merge.annot2, aes(x = stype, y = log2(chromocounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~plate, nrow = 1)
  
ggplot(dat.merge.annot2, aes(x = stype, y = log2(chromocounts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  facet_wrap(~plate, nrow = 1)



dev.off()