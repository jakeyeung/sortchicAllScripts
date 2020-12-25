# Jake Yeung
# Date of Creation: 2020-12-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/10-check_glmpca_downstream.othermarks.use_same_annot_file.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(JFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load glmpca  -------------------------------------------------------------

# jmark <- "H3K4me1"

jalpha <- 1

nitervec <- c("500", "1000")
binskeepvec <- c("0", "500", "1000")
# niter <- "500"
# binskeep <- 1000
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K27me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"
indir.glmpca <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/glmpca_outputs_BM_AllMerged3_plate.same_annot_file_rerun"
dir.create(outdir)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

jdate <- "2020-12-03"
for (niter in nitervec){
  for (binskeep in binskeepvec){
    
    jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
    
    # indir.glmpca <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs"
    
    for (jmark in jmarks){
      
      outbase <- file.path(outdir, paste0("glmpca_plate_effects_output.", jmark, ".", jsuffix, ".", jdate))
      outpdf <- paste0(outbase, ".pdf")
      outtxt.glmpca <- paste0(outbase, ".glmpca.txt")
      outtxt.lda <- paste0(outbase, ".lda.txt")
      
      fname <- paste0("glmpca.", jmark, ".", jsuffix, ".RData")
      inf.glmpca <- file.path(indir.glmpca, fname)
      
      if (file.exists(outpdf)){
         print(paste("Skipping writing of", outpdf))
         next
      }
      
      if (!file.exists(inf.glmpca)){
        print(paste("File does not exist, skipping:", inf.glmpca))
        next
      }
      
      pdf(outpdf, useDingbats = FALSE)
      
      print(paste("Opening:", inf.glmpca))
      
      load(inf.glmpca, v=T)
      
      dat.glmpca <- DoUmapAndLouvain(glm.out$factors, jsettings)
      
      
      ggplot(dat.glmpca, mapping = aes(x = umap1, y = umap2, color = louvain)) + 
        geom_point() + 
        scale_color_manual(values = cbPalette) + 
        theme_bw() + 
        ggtitle(jmark, "GLMPCA") + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      
      # Compare with LDA ?  -----------------------------------------------------
      
      # indir.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins"
      indir.lda <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file")
      # fname.lda <- paste0("lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj")
      fname.lda <- paste0("lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.Robj")
      inf.lda <- file.path(indir.lda, fname.lda)
      
      load(inf.lda, v=T)
      
      tm.result <- posterior(out.lda)
      dat.lda <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)
      
      ggplot(dat.lda, mapping = aes(x = umap1, y = umap2, color = louvain)) + 
        geom_point() + 
        scale_color_manual(values = cbPalette) + 
        theme_bw() + 
        ggtitle(jmark, "LDA") + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      
      # Add celltype info -------------------------------------------------------
      
      # indir.annot <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04"
      # fname.annot <- paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt")
      # dat.annot <- fread(file.path(indir.annot, fname.annot)) %>%
      #   rowwise() %>%
      #   mutate(plate = ClipLast(cell, jsep = "_"),
      #          experi = ClipLast(plate, jsep = "-"))
      
      inf.annot <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.", jmark, ".2020-11-18.dupfilt.txt"))
      dat.annot <- fread(inf.annot)
      
      dat.glmpca.annot <- left_join(dat.glmpca, subset(dat.annot, select = c(cell, cluster, plate)))
      dat.lda.annot <- left_join(dat.lda, subset(dat.annot, select = c(cell, cluster, plate)))
      
      m.glmpca <- ggplot(dat.glmpca.annot, aes(x = umap1, y = umap2, color = cluster)) + 
        geom_point(alpha = jalpha) +  
        scale_color_manual(values = cbPalette) + 
        ggtitle(jmark, "GLMPCA") + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      
      m.glmpca.plates <- ggplot(dat.glmpca.annot, aes(x = umap1, y = umap2, color = cluster)) + 
        geom_point(alpha = jalpha) +  
        scale_color_manual(values = cbPalette) + 
        ggtitle(jmark, "GLMPCA") + 
        facet_wrap(~plate) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      
      m.lda <- ggplot(dat.lda.annot, aes(x = umap1, y = umap2, color = cluster)) + 
        geom_point(alpha = jalpha) +  
        scale_color_manual(values = cbPalette) + 
        ggtitle(jmark, "LDA") + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 
      
      m.lda.plates <- ggplot(dat.lda.annot, aes(x = umap1, y = umap2, color = cluster)) + 
        geom_point(alpha = jalpha) +  
        scale_color_manual(values = cbPalette) + 
        ggtitle(jmark, "LDA") + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
        facet_wrap(~plate)
      
      
      multiplot(m.glmpca, m.lda, cols = 2)
      print(m.glmpca.plates)
      print(m.lda.plates)
      # saveoutput
      fwrite(dat.glmpca.annot, file = outtxt.glmpca, quote = FALSE, sep = "\t", col.names = TRUE, na = "NA")
      fwrite(dat.lda.annot, file = outtxt.lda, quote = FALSE, sep = "\t", col.names = TRUE, na = "NA")
      dev.off()
    }
    
  }
}




# Do stuff  ---------------------------------------------------------------






