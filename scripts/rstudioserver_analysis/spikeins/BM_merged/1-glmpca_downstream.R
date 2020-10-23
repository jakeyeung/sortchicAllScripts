# Jake Yeung
# Date of Creation: 2020-10-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/1-LDA_downstream.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)


library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

SplitGetLast <- function(x, jsplit = "-"){
  # get last element after splitting
  xsplit <- strsplit(x, jsplit)[[1]]
  xnew <- xsplit[[length(xsplit)]]
  return(xnew)
}

GetRepBM <- function(experiname){
  if (grepl("rep3", experiname)){
    return("rep3")
  } else if (grepl("rep2", experiname)){
    return("rep2")
  } else {
    print("Neither rep2 nor rep3 in name: Returning rep2")
    return("rep2")
  }
}


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115")


# Load data  --------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein_BMround2all"
assertthat::assert_that(dir.exists(indir))

inmain.lda <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all")
indir.lda <- file.path(inmain.lda)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

jmaxiter <- 1000
jtol <- as.character("1e-6")

# jmaxiter <- 1000
# jtol <- as.character("1e-6")

jsuffix <- ".filt_0.15_0.95_counts_and_l2r"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all/glmpca_outputs"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("glmpca_output.maxiter_", jmaxiter, "tol_", jtol, ".pdf"))



pdf(outpdf, useDingbats = FALSE)

dat.umap.lst <- lapply(jmarks, function(jmark){
  
  print(jmark)
  inf <- file.path(indir, paste0("count_mat.", jmark, jsuffix, ".glmpcaout.penalty_1.maxiter_", jmaxiter, ".stochastic.avagrad.tol_", jtol, ".devfilt_5000.by_plate.RData"))
  print(basename(inf))
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  dat.umap.long <- DoUmapAndLouvain(glmpcaout$factors, jsettings = jsettings)
  
  dat.umap.long$dim1 <- glmpcaout$factors[, 1]
  dat.umap.long$dim2 <- glmpcaout$factors[, 2]
  dat.umap.long$dim3 <- glmpcaout$factors[, 3]
  dat.umap.long$dim4 <- glmpcaout$factors[, 4]
  
  dat.umap.long$mark <- jmark
  
  return(dat.umap.long)
})


# jmark <- jmarks[[1]]
for (jmark in jmarks){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}



dat.umap.long <- dat.umap.lst %>% bind_rows() %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"),
         plate = as.numeric(SplitGetLast(experi, jsplit = "-")), 
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))

for (jmark in jmarks){
  m <- ggplot(dat.umap.long %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = as.character(plate))) + 
    geom_point(alpha = 0.33)  +  
    ggtitle(jmark) + 
    facet_wrap(~stype) + 
    scale_color_manual(values = cbPalette) +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.umap.long %>% filter(mark == jmark), aes(x = dim1, y = dim2, color = stype)) + 
    geom_point(alpha = 0.33)  +  
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


# Check chromo to spikeints -----------------------------------------------

# inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt.VAN5046/spikein_info.txt"
inf.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/spikein_info_BM_round2_all.txt"
dat.spikeins.mat <- fread(inf.spikeins)

dat.merge <- left_join(dat.umap.long, subset(dat.spikeins.mat, select = c(samp, spikeincounts, chromocounts)), by = c("cell" = "samp"))

for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = stype)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    facet_wrap(~stype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = dim1, y = dim2, color = log2(chromocounts / spikeincounts))) + 
    geom_point(alpha = 0.33)  +  
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    facet_wrap(~stype) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}

for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = dim1, y = dim2, color = log2(chromocounts / spikeincounts))) + 
    geom_point(alpha = 0.33)  +  
    ggtitle(jmark) + 
    scale_color_viridis_c() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
}


# Check variances ---------------------------------------------------------

# load LDAs calculate variances


dat.vars.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  outname <- file.path(paste0("lda_outputs.count_mat.", jmark, jsuffix, ".K-30.binarize.FALSE"), paste0("ldaOut.count_mat.", jmark, jsuffix, ".K-30.Robj"))
  inf <- file.path(indir.lda, outname)
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  tm.result <- posterior(out.lda)
  dat.impute.log <- t(log2(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  return(dat.var)
})

dat.merge.var <- left_join(dat.merge, dat.vars.lst %>% bind_rows(), by = "cell")

for (jmark in jmarks){
  m <- ggplot(dat.merge.var %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_viridis_c(direction = -1, na.value = "grey85") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


dev.off()






