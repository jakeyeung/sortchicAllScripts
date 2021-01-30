# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/5-write_cell_cluster_tables_updated_umap.R
# Write cell cluster tables, H3K27me3 stays the same because there is no old cells

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
stypecols <- c("grey", "red", "blue")

ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")

jsize <- 0.6
jalpha <- 0.8

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28"
dir.create(outdir)

# Load metas --------------------------------------------------------------


indir.metas <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned")
assertthat::assert_that(dir.exists(indir.metas))

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  fread(file.path(indir.metas, fname))
})


# Load glmpca -------------------------------------------------------------


# jmark <- "H3K9me3"
niter <- "100"
inf.glmpca.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/BM_from_matadj/glmpca.", jmark, ".from_matadj.platename_jrep.szname_none.niter_", niter, ".RData"))
    assertthat::assert_that(file.exists(inf.glmpca))
  } else {
    inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_plate.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData"))
  }
  return(inf.glmpca)
})

glmpca.factors <- lapply(inf.glmpca.lst, function(inf){
  load(inf, v=T)
  return(glm.out$factors)
})


# Update umaps  -----------------------------------------------------------

dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  dat.umap <- DoUmapAndLouvain(glmpca.factors[[jmark]], jsettings = jsettings) 
  dat.umap.annot <- left_join(dat.umap, dat.metas[[jmark]] %>% dplyr::select(c(-umap1, -umap2, -louvain)))
  return(dat.umap.annot)
})


# Flip x axes for K4me3 K9me3  --------------------------------------------

dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  jtmp <- dat.umap.annot.lst[[jmark]]
  if (jmark %in% c("H3K4me3", "H3K9me3")){
    jtmp$umap1 <- -1 * jtmp$umap1
  }
  return(jtmp)
})




# Write outputs -----------------------------------------------------------

# glmpca factors and umaps 
for (jmark in jmarks){
  print(jmark)
  # write glmpca factors
  outfactors <- file.path(outdir, paste0("glmpca_factors_batch_corrected.", jmark, ".", Sys.Date(), ".txt"))
  outannots<- file.path(outdir, paste0("metadata_batch_corrected.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(glmpca.factors[[jmark]], file = outfactors, quote = FALSE, sep = "\t", row.names = TRUE)
  fwrite(dat.umap.annot.lst[[jmark]], file = outannots, quote = FALSE, sep = "\t")
}


# Write pdfs  -------------------------------------------------------------

outpdf <- file.path(outdir, paste0("plots_batch_corrected.", jmark, ".", Sys.Date(), ".pdf"))
pdf(outpdf, useDingbats = FALSE)

# show no batch efects
for (jmark in jmarks){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(size = jsize, alpha = jalpha) + 
    theme_minimal(4) + 
    ggtitle(jmark) + 
    facet_wrap(~jrep) + 
    scale_color_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

# color by cluster
for (jmark in jmarks){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(size = jsize, alpha = jalpha) + 
    theme_minimal(4) + 
    ggtitle(jmark) + 
    scale_color_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

# color by batch
for (jmark in jmarks){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = stypecol)) + 
    geom_point(size = jsize, alpha = jalpha) + 
    theme_minimal(4) + 
    ggtitle(jmark) + 
    scale_color_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


# show legend: cluster
for (jmark in jmarks){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point(size = jsize, alpha = jalpha) + 
    theme_minimal(4) + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

# show legend: batch
for (jmark in jmarks){
  print(jmark)
  m <- ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = batch)) + 
    geom_point(size = jsize, alpha = jalpha) + 
    theme_minimal(4) + 
    ggtitle(jmark) + 
    scale_color_manual(values = stypecols, na.value = "grey85") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


dev.off()