# Jake Yeung
# Date of Creation: 2020-12-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3b-make_umaps_with_TSS_cuts.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(DescTools)
library(ggrastr)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jmarks.act <- c("H3K4me1", "H3K4me3"); names(jmarks.act) <- jmarks.act

pdfout <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/famous_genes_on_umap_raw_cuts.ListFromPeter.pdf"

pdf(pdfout, useDingbats = FALSE)

# Load UMAPs --------------------------------------------------------------

indir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final")

dat.metas <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt"))
  dat.meta <- fread(inf)
})


# Load raw cuts  ----------------------------------------------------------

inf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/raw_cuts_TSS_three_marks.rds"
dat.raw.lst <- readRDS(inf.tss)
dat.raw.lst <- lapply(dat.raw.lst, function(jdat){
  sweep(x = jdat, MARGIN = 2, STATS = colSums(jdat), FUN = "/")
})

rnames.act <- lapply(jmarks.act, function(jmark){
  rnames <- rownames(dat.raw.lst[[jmark]])
}) %>%
  unlist() %>%
  unique()  %>%
  gtools::mixedsort()

coords.act <- sapply(rnames.act, function(x) paste("chr", strsplit(x, split = ";")[[1]][[1]], sep = ""))
coord2rname <- hash::hash(coords.act, rnames.act)

rownames(dat.raw.lst$H3K27me3) <- sapply(rownames(dat.raw.lst$H3K27me3), function(x) AssignHash(x = x, jhash = coord2rname, null.fill = x))



#  Get famous genes  -------------------------------------------------------

# famous.genes.lst <- 
#   list("Eryths" = c("Tal1", "Aqp1", "Gata1", "Comt", "Sphk1", "Hbb-bs", "Hbb-bt", "Sox6"),
#        "Bcells" = c("Pax5", "Ebf1", "Cd79a", "Cd79b", "Snn", "Blnk", "Cd72", "Blk", "Kzf3", "Cd19"), 
#        "NKs" = c("Stat4", "Tcf7", "Tbx21", "Cd3d", "Gimap4"), 
#        "Granulocytes" = c("Cxcr2", "Mmp9", "Cd177", "Ncf4", "S100a9", "S100a8", "Ncf1", "Cebpe"),
#        "Basophils" = c("Il4", "Il1r1", "Arap3", "Cpa3"), 
#        "pDCs" = c("Irf8", "Selpg", "Itgb7", "Ccr9", "Unc93b1"), 
#        "DCs" = c("Ccl9", "Apoe", "Nlrp3", "Csf1r"), 
#        "HSPCs" = c("Hoxa9", "Hoxa7", "Hoxa3", "Meis1", "Runx2", "Kit", "Hlf", "Hoxa10", "Erg", "Cd34", "Hoxa6", "Gata3", "Hoxb4"))

famous.genes.lst <-
  list("Eryths" = c("Hbb-bs"),
       "Bcells" = c("Vpreb1", "Vpreb2", "Cd79a"),
       "NKs" = c("Gzma"),
       "Granulocytes" = c("Itgb2l", "Elane", "Fcnb", "Mpo", "S100a9", "S100a8", "Prtn3", "Lcn2", "Lrg1"),
       "Basophils" = c("Prg3", "Epx", "Prg2"),
       "pDCs" = c("Siglech"),
       "DCs" = c("Cst3", "Xcr1"),
       "HSPCs" = c("Hoxa9", "Hoxa3", "Meis1", "Kit", "Hlf", "Hoxa10", "Hoxa6"))


jgene <- "Tal1"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
lapply(jmarks, function(jmark){
  m <- ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) + 
    geom_point(alpha = 0.2) + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_identity(guide = "legend", 
                         labels = unique(dat.metas[[jmark]]$cluster)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
})


for (jgene in unlist(famous.genes.lst)){
  print(jgene)
  dat.tss.lst <- lapply(jmarks, function(jmark){
    print(jmark)
    rname <- grep(jgene, rownames(dat.raw.lst[[jmark]]), value = TRUE)
    if (length(rname) == 0){
      print("Skipping, nothing found")
      print(rname)
      return(NULL)
    } else {
      rname <- rname[[1]]
    }
    # print(rname)
    dat.tss <- data.frame(cuts = dat.raw.lst[[jmark]][rname, ], cell = colnames(dat.raw.lst[[jmark]]), stringsAsFactors = FALSE) %>%
      left_join(., dat.metas[[jmark]]) %>%
      arrange(cuts) %>%
      ungroup() %>%
      mutate(cuts.win = Winsorize(log2(cuts + 1), probs = c(0, 0.99)),
             gene = jgene,
             tssname = rname)
    return(dat.tss)
  })
  if (is.null(dat.tss.lst[[1]])){
    print("Skipping, nothing found2")
    next
  }
  
  mlst <- lapply(jmarks, function(jmark){
    jrname <- dat.tss.lst[[jmark]]$tssname[[1]]
    m <- ggplot(dat.tss.lst[[jmark]], aes(x = umap1, y = umap2, color = cuts.win)) + 
      geom_point_rast() + 
      scale_color_viridis_c() + 
      theme_minimal() + 
      ggtitle(paste(jmark, jgene, jrname)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  })
  print(mlst)
  # JFuncs::multiplot(mlst$H3K4me3, mlst$H3K4me1, mlst$H3K27me3, cols = 3)
}

dev.off()