# Jake Yeung
# Date of Creation: 2021-01-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/8-MARA_downstream_k9dynamicbins.50kbbins.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)
library(ggrepel)


hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jsize <- 1.5
zscores.cutoff <- 1.25

# Get emtas ---------------------------------------------------------------


ctypes <- c("Eryths", "Bcells", "NKs", "pDCs", "Granulocytes", "DCs", "Basophils", "HSPCs")
# indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned")
indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")
assertthat::assert_that(dir.exists(indir.meta))
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  # fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  fname.tmp <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.meta, fname.tmp)
  print(inf)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})


# jmark <- "H3K4me1"
for (jmark in jmarks){
  
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/MARA_output_BM_k9dynamicbins"
  fname <- paste0("MARA_output_BM_k9dynamicbinsfilt.", jmark, ".50kbins.zscore_", zscores.cutoff, ".dotsize_", jsize, ".", Sys.Date(), ".pdf")
  outpdf <- file.path(outdir, fname)
  
  pdf(outpdf, useDingbats = FALSE)
  
  
  
  
  # Load MARA output --------------------------------------------------------
  
  # jsuffix2 <- paste0("BM_", jmark, ".BM_AllMerged3.glmpca_plate.bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt.cleanuprows.same_annot_file")
  # jsuffix2 <- paste0("ldaOut.count_mat_by_peaks.", jmark, ".overlap_dynamic_k9_bins.K-30.keepNbins_0-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/ldaOut.count_mat_by_peaks.H3K4me1.overlap_dynamic_k9_bins.K-30.keepNbins_0")
  jsuffix2 <- paste0("ldaOut.count_mat_k9_dynamic_bins_50kb.", jmark, ".2021-01-28.K-30.keepNbins_0-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.H3K9me3/ldaOut.count_mat_k9_dynamic_bins_50kb.", jmark, ".2021-01-28.K-30.keepNbins_0")
  
  # indir.mara <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/", jmark, "/mara_output.lessmemory/", jsuffix2, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/", jsuffix2))
  indir.mara <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/", jmark, "/mara_output/", jsuffix2))
  print(indir.mara)
  assertthat::assert_that(dir.exists(indir.mara))
  
  mara.out <- LoadMARA(indir.mara, make.cnames = FALSE)
  
  act.mat.clean <- t(subset(mara.out$act.mat, select = -motif))
  colnames(act.mat.clean) <- mara.out$act.mat$motif
  act.mat.clean.dat <- data.frame(cell = rownames(act.mat.clean), act.mat.clean, stringsAsFactors = FALSE) %>% 
    ungroup() %>%
    mutate(cell = gsub("\\.", "-", cell))
  
  
  motifs.keep <- subset(mara.out$zscores, zscore > zscores.cutoff)$motif
  
  m.zscores <- mara.out$zscores %>%
    ungroup() %>% 
    mutate(rnk = seq(length(motif)),
           motiflab = ifelse(zscore > zscores.cutoff * 1.5, toupper(motif), NA)) %>%
    ggplot(., aes(x = rnk, y = zscore, label = motiflab)) + 
    geom_point() + 
    geom_text_repel() + 
    theme_bw() + 
    geom_hline(yintercept = zscores.cutoff, linetype = "dotted") + 
    xlab("Rank") + 
    ylab("Motif zscore") + 
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.zscores)
  
  m.zscores2 <- mara.out$zscores %>%
    ungroup() %>% 
    mutate(rnk = seq(length(motif))) %>%
    ggplot(., aes(x = zscore)) + 
    geom_density(fill = "red", alpha = 0.25) + 
    geom_vline(xintercept = zscores.cutoff, linetype = "dotted") + 
    theme_bw() + 
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.zscores2)
  
  
  
  # add motifs to UMAP  -----------------------------------------------------
  
  
  
  dat.merge.motifs <- left_join(dat.metas[[jmark]], act.mat.clean.dat, by = "cell")
  
  print(head(motifs.keep))
  
  jmotif <- "Zfp110"
  jmotif <- "Zfp110"
  jmotif <- "Zfp711"
  jmotif <- "Spic"
  jmotif <- "Gfi1"
  
  for (jmotif in motifs.keep){
    jzscore <- signif(subset(mara.out$zscores, motif == jmotif)$zscore, digits = 2)
    (jtitle <- paste(jmotif, "Zscore:", jzscore))
    
    # pdf("/home/jyeung/hub_oudenaarden/jyeung/tmp/motiftest.pdf", useDingbats = FALSE)
    m <- ggplot(dat.merge.motifs, aes_string(x = "umap1", y = "umap2", color = jmotif)) + 
      # geom_point(size = 0.75) + 
      geom_point(size = jsize) + 
      ggtitle(jtitle) + 
      theme_bw() + 
      scale_color_viridis_c() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # m <- PlotXYWithColor(dat.merge.motifs, xvar = "umap1", yvar = "umap2", cname = jmotif, jtitle = jtitle, cont.color = TRUE, jsize = 0.5) + scale_color_viridis_c()
    print(m)
  }
  
  
  
  # Add heatmap -------------------------------------------------------------
  
  
  jmat <- as.matrix(subset(act.mat.clean.dat, select = -cell))
  rownames(jmat) <- act.mat.clean.dat$cell
  
  # jmat[dat.metas$H3K9me3$cell, ][1:5, 1:5]
  
  cells.ordered <- dat.metas[[jmark]]$cell
  jsub <- jmat[cells.ordered, motifs.keep]
  
  
  library(heatmap3)
  
  colvec <- dat.metas[[jmark]]$clustercol
  
  jmeth <- "ward.D2"
  
  par(mfrow=c(1,1), mar=c(1,1,1,1), mgp=c(3, 1, 0), las=0)
  hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                     main = jmark, 
                     # ColSideColors = rep("blue", ncol(jsub)), 
                     # ColSideColors = FALSE,
                     RowSideColors = colvec, 
                     # RowSideColors = rep("red", nrow(jsub)), 
                     RowSideLabs = "celltype", 
                     labRow = FALSE, scale = "column", revC = TRUE,
                     distfun = dist, hclustfun = hclust, method = jmeth)
  
  hm.out <- heatmap3(jsub, margins = c(5, 8), cexCol = 0.35, Colv = TRUE, Rowv = NA, 
                     main = jmark, 
                     # ColSideColors = rep("blue", ncol(jsub)), 
                     # ColSideColors = FALSE,
                     RowSideColors = colvec, 
                     # RowSideColors = rep("red", nrow(jsub)), 
                     RowSideLabs = "celltype", 
                     labRow = FALSE, scale = "column", revC = FALSE,
                     distfun = dist, hclustfun = hclust, method = jmeth)
  
  
  hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                     main = jmark, 
                               # ColSideColors = rep("blue", ncol(jsub)), 
                               # ColSideColors = FALSE,
                               ColSideColors = colvec, 
                               # RowSideColors = rep("red", nrow(jsub)), 
                               ColSideLabs = "celltype", 
                               labCol = FALSE, scale = "row", revC = FALSE, 
                               distfun = dist, hclustfun = hclust, method = jmeth)
  
  hm.out.transpose <- heatmap3(t(jsub), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE, 
                     main = jmark, 
                               # ColSideColors = rep("blue", ncol(jsub)), 
                               # ColSideColors = FALSE,
                               ColSideColors = colvec, 
                               # RowSideColors = rep("red", nrow(jsub)), 
                               ColSideLabs = "celltype", 
                               labCol = FALSE, scale = "row", revC = TRUE,
                               distfun = dist, hclustfun = hclust, method = jmeth)
  
  
  
  
  
  dev.off()
}




