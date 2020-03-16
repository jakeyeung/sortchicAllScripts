# Jake Yeung
# Date of Creation: 2020-03-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/1-explore_MARA_outputs_marks_merged.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)
library(JFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load LDA and glm --------------------------------------------------------

# jmark <- "H3K4me3"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3")
names(jmarks) <- jmarks

# jmdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_output/NoIntercept_ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged/ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30")
jmdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_output/NoIntercept_ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged/ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30")
jmdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_output/AddIntercept_ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.marks_merged/ldaOut.marks_merged.BM_AllMerged.merged_by_clusters_no_NAs.K-30")
assertthat::assert_that(dir.exists(jmdir))
dat.umap.long.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  
  inf.glmpca <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData")
  load(inf.glmpca, v=T)
  
  
  # Plot outputs ------------------------------------------------------------
  
  dat.umap.long <- DoUmapAndLouvain(glm.out$factors, jsettings)
  
  dat.umap.long <- dat.umap.long %>%
    rowwise() %>%
    mutate(experi = ClipLast(x = as.character(cell), jsep = "_"),
           cond = GetCondFromSamp(as.character(cell), mark = jmark)) %>%
    ungroup() %>%
    mutate(cond = factor(x = cond, levels = c("Unenriched", "Linneg", "StemCell")))
  return(dat.umap.long)
})


mara.out <- LoadMARA(mdir = jmdir, make.cnames = FALSE)

act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
colnames(act.mat.clean) <- mara.out$act.mat$motif
act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
  ungroup() %>%
  mutate(cell = gsub("\\.", "-", cell))

# Merge datumaps to mara  -------------------------------------------------

names(dat.umap.long.lst) <- jmarks

# split matrix by mark and then plot

mat.split <- lapply(jmarks, function(jmark){
  subset(act.mat.clean.dat, grepl(jmark, cell))
})

dat.merge <- lapply(jmarks, function(jmark){
  left_join(dat.umap.long.lst[[jmark]], mat.split[[jmark]], by = "cell")
})

# Plot output -------------------------------------------------------------

jmark <- "H3K4me1"
jmotif <- "Gata1"
jzscore <- signif(subset(mara.out$zscores, motif == jmotif)$zscore, digits = 2)
jtitle <- paste(jmotif, "Zscore:", jzscore)

m.lst <- lapply(jmarks, function(jmark){
  m <- PlotXYWithColor(dat.merge[[jmark]], xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE) + scale_color_viridis_c()
  return(m)
})

multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)





