# Jake Yeung
# Date of Creation: 2020-03-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/1-explore_MARA_outputs_centered_vs_uncentered.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load LDA and glm --------------------------------------------------------

# jmark <- "H3K4me3"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3")

for (jmark in jmarks){
  print(jmark)
  ntopics <- "30"
  jbinskeep <- 250
  jsubdir <- "centered"
  
  inf.glmpca <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData")
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains/lda_outputs.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".binarize.FALSE/ldaOut.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".Robj")
  inf.lda.bins <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
  
  load(inf.glmpca, v=T)
  # load(inf.lda, v=T)
  
  outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/mara_outputs/outputs_uncentered"
  outname <- paste0(ClipLast(basename(inf.glmpca), jsep = "\\.", jsep.out = "."), ".pdf")
  outpdf <- file.path(outdir, outname)
  
  # Load MARA  --------------------------------------------------------------
  
  # if (jsubdir == "centered"){
  #   # jmdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/NoIntercept_countmat_PZ_fromHiddenDomains_H3K4me1.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.H3K4me1/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
  # } else if (jsubdir == "uncentered"){
  #   jmdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/AddIntercept_countmat_PZ_fromHiddenDomains_H3K4me1.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.H3K4me1/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
  # } else {
  #   warning("Msut be jsubdir centered or uncentered found:", jsubdir)
  # }
  
  jmdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/AddIntercept_countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.", jmark, "/countmat_PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_250.ntopics_30.2020-02-11")
  print(jmdir)
  assertthat::assert_that(dir.exists(jmdir))
  
  mara.out <- LoadMARA(mdir = jmdir, make.cnames = FALSE)
  
  print(head(mara.out$zscores))
  
  
  # Plot outputs ------------------------------------------------------------
  
  
  dat.umap.long <- DoUmapAndLouvain(glm.out$factors, jsettings)
  
  dat.umap.long <- dat.umap.long %>%
    rowwise() %>%
    mutate(experi = ClipLast(x = as.character(cell), jsep = "_"),
           cond = GetCondFromSamp(as.character(cell), mark = jmark)) %>%
    ungroup() %>%
    mutate(cond = factor(x = cond, levels = c("Unenriched", "Linneg", "StemCell")))
  
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cond)) + geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
  
  
  # Plot activities ---------------------------------------------------------
  
  mara.out$act.mat[1:5, 1:5]
  
  act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
  colnames(act.mat.clean) <- mara.out$act.mat$motif
  act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
    ungroup() %>%
    mutate(cell = gsub("\\.", "-", cell))
  
  dat.merge <- left_join(dat.umap.long, act.mat.clean.dat, by = "cell")
  
  # jmotif <- "Tal1"
  # PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jmotif, cont.color = TRUE) + scale_color_viridis_c()
  # jmotif <- "Ebf1"
  # jmotif <- "Hoxb5"
  # PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jmotif, cont.color = TRUE) + scale_color_viridis_c()
  # plot top zscores
  
  zscores.sub <- mara.out$zscores %>%
    ungroup() %>%
    filter(motif != "intercept") %>%
    top_n(x = ., n = 50, wt = zscore)
  
  # plot the top 
  
  jmotifs <- zscores.sub$motif
  
  
  pdf(file = outpdf, useDingbats = FALSE)
  for (jmotif in jmotifs){
    jzscore <- subset(mara.out$zscores, motif == jmotif)$zscore
    jtitle <- paste(jmotif, "Zscore:", signif(jzscore, digits = 2))
    m <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE) + scale_color_viridis_c()
    print(m)
    m.split <- m + facet_wrap(~cond)
  }
  dev.off()
  
}



