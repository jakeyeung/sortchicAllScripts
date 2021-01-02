# Jake Yeung
# Date of Creation: 2020-12-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/11-tweak_UMAP_parameters_write_glmpca_factors.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(hash)
library(igraph)
library(umap)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load UMAPs blow them up  ------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

dat.metas <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt"))
  fread(inf.meta)
})

library(ggforce)
mlst <- lapply(jmarks, function(jmark){
  ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(size = 0.2, position = position_jitternormal(sd_x = 0.5, sd_y = 0.5), alpha = 1) + 
    ggtitle(jmark) + 
    theme_bw() + 
    scale_color_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)

# Try loading from glmpca  ------------------------------------------------


niter <- "500"
niter2 <- "500"
# niter2 <- "1000"
binskeep <- "0"
binskeep2 <- "1000"
jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
jsuffix2 <- paste0("bincutoff_0.binskeep_", binskeep2, ".byplate.szname_none.niter_", niter2, ".reorder_rownames.dupfilt")

infs.lst <- lapply(jmarks, function(jmark){
  if (jmark == "H3K27me3"){
    # inf <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep1rep2rep3reseq.peaks.varfilt/glmpca.H3K27me3.bincutoff_0.binskeep_0.platename_jrep.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData")
    inf <- file.path(hubprefix, "jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep2rep3reseq.peaks.varfilt/glmpca.H3K27me3.bincutoff_0.binskeep_0.platename_plate.szname_none.niter_500.reorder_rownames.dupfilt.suffix_peaks.RData")
  # } else if (jmark == "H3K9me3"){
  } else if (jmark == c("H3K9me3")){
    print(paste(jmark, "special"))
    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix2, ".RData"))
  } else {
    # use same annot file rerun because basophils are missing for K4me1 otherwise
    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun/glmpca.", jmark, ".", jsuffix, ".RData"))
  }
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

glmout.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- infs.lst[[jmark]]
  load(inf, v=T)
  return(glm.out$factors)
})

# Write factors to output -------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final_glmpca_tweak_umap.H3K27me3_rep2rep3reseq"
dir.create(outdir)

for (jmark in jmarks){
  outftmp <- file.path(outdir, paste0("glmpca_factors.cleaned.", jmark, ".", Sys.Date(), ".txt"))
  # write.table(glmout.lst[[jmark]], file = outftmp, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  fwrite(glmout.lst[[jmark]], file = outftmp, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}


# Vary spread and mindist -------------------------------------------------

# vary spreads
outpdf <- file.path(outdir, paste0("umaps_tweaks.H3K27me3_rep2rep3reseq", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

spreads <- c(5, 8)
mindists <- c(0.1)

# vary spreads
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

for (spread in spreads){
  print('spread')
  print(spread)
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$spread <- spread
  jsettings$random_state <- 123

  dat.umap.lst <- mclapply(jmarks, function(jmark){
  # dat.umap.lst <- lapply(jmarks, function(jmark){
    print(jmark)
    glmout <- glmout.lst[[jmark]]
    dat.umap <- DoUmapAndLouvain(glmout, jsettings = jsettings)
  }, mc.cores = 4)
  # })
  
  mlst <- lapply(jmarks, function(jmark){
    m <- ggplot(dat.umap.lst[[jmark]] %>% left_join(., dat.metas[[jmark]] %>% dplyr::select(c("cell", "cluster"))), 
           aes(x = umap1, y = umap2, color = cluster)) + 
      geom_point(size = 0.2, alpha = 0.2) + 
      ggtitle(jmark, paste(paste0("spread=", jsettings$spread), paste0("\nmindist=", jsettings$min_dist), sep = ",")) + 
      theme_bw() + 
      scale_color_manual(values = cbPalette) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
    # test 
    mcheck <- ggplot(dat.umap.lst[[jmark]] %>% left_join(., dat.metas[[jmark]] %>% dplyr::select(c("cell", "cluster", "jrep"))), 
                     aes(x = umap1, y = umap2, color = cluster)) + 
      geom_point(size = 0.8, alpha = 0.8) + 
      facet_wrap(~jrep) + 
      ggtitle(jmark, paste(paste0("spread=", jsettings$spread), paste0("\nmindist=", jsettings$min_dist), sep = ",")) + 
      theme_bw() + 
      scale_color_manual(values = cbPalette) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(mcheck)
    
    return(m)
  })
  JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)
}


for (mindist in mindists){
  print('mindist')
  print(mindist)
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- mindist
  jsettings$spread <- 1
  jsettings$random_state <- 123
  
  dat.umap.lst <- mclapply(jmarks, function(jmark){
    glmout <- glmout.lst[[jmark]]
    dat.umap <- DoUmapAndLouvain(glmout, jsettings = jsettings)
  }, mc.cores = 4)
  
  mlst <- lapply(jmarks, function(jmark){
    ggplot(dat.umap.lst[[jmark]] %>% left_join(., dat.metas[[jmark]] %>% dplyr::select(c("cell", "cluster"))), 
           aes(x = umap1, y = umap2, color = cluster)) + 
      geom_point(size = 0.2, alpha = 0.2) + 
      ggtitle(jmark, paste(paste0("spread=", jsettings$spread), paste0("\nmindist=", jsettings$min_dist), sep = ",")) + 
      theme_bw() + 
      scale_color_manual(values = cbPalette) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  })
  JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)
  
}

dev.off()


# Write output so people can play with it  --------------------------------



