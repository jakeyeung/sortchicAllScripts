# Jake Yeung
# Date of Creation: 2022-04-09
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8-downstream_LL_celltyping_repressive_cleaned_update_ctypes_from_LDA_newannot_eucdist.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k4me1", "k9me3", "k27me3", "k4me3"); names(jmarks) <- jmarks
jsuffix <- "dynamicbins"
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_newannot"
outdir <- indir

# jmark <- jmarks[[1]]


# Get colors --------------------------------------------------------------

dat.meta.colors <- fread("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap/metadata.k4me1.2022-03-01.txt") %>%
  dplyr::select(ctype.from.LL, colcode) 
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

# replace monocyt ecolor with #d69941
colors.hash <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)
colors.hash[["Monocytes"]] <- "#d69941"
colors.hash[["CLP"]] <- "#ACD4E5"


for (jmark in jmarks){
  print(jmark)
  
  outmeta <- file.path(outdir, paste0("metadata_reannotate_from_LLmat_fix_ctypes_newannot_eucdist_", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  outpdf <- file.path(outdir, paste0("plots_reannotate_from_LLmat_fix_ctypes_newannot_eucdist_", jsuffix, ".", jmark, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  inf <- file.path(indir, paste0("LLmat_by_batch_dynamicbins_", jmark, ".2022-04-10.RData"))
  
  load(inf, v=T)
  
  if (jmark %in% c("k4me1", "k4me3")){
    indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap"
    inf.meta <- file.path(indir.meta, paste0("metadata.", jmark, ".2022-03-01.txt"))
    dat.meta.final <- fread(inf.meta) %>%
      dplyr::select(-ctype.from.LL)
  } else if (jmark %in% c("k27me3", "k9me3")) {
    inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_check_eryths/metadata_celltyping_", jmark, ".dynamicbins.2022-04-06.txt")
    dat.meta.final <- fread(inf.meta)
  }
  
  inf.meta.newannot <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata_new_prog_annots/metadata_newProgenitors.", jmark, ".2022-03-01.txt")
  dat.meta.newannot <- fread(inf.meta.newannot) %>%
    dplyr::select(cell, ctype)
  
  dat.meta.newannot.notnan <- subset(dat.meta.newannot, ctype != "NaN") 
  
  reannotate.ground.truth <- hash::hash(dat.meta.newannot.notnan$cell, dat.meta.newannot.notnan$ctype)
  
  # edit dat.meta.final already
  dat.meta.final <- dat.meta.final %>%
    rowwise() %>%
    mutate(ctype = AssignHash(x = cell, jhash = reannotate.ground.truth, null.fill = ctype))
  
  closest.cells.best <- c(closest.cells.best.old, closest.cells.best.new)
  closest.cells.hash <- hash::hash(names(closest.cells.best), closest.cells.best)
  
  # dat.meta.newannot.reannotate <- dat.meta.newannot %>%
  #   rowwise() %>%
  #   mutate(ctype.from.LL = AssignHash(x = cell, jhash = closest.cells.hash, null.fill = ctype),
  #          ctype.from.LL = ifelse(ctype.from.LL == "CMP/GMP", "CMP", ctype.from.LL),
  #          ctype.from.LL = ifelse(ctype.from.LL == "LT", "HSCs", ctype.from.LL), 
  #          ctype.from.LL = ifelse(ctype.from.LL == "ST", "HSCs", ctype.from.LL))
  #          # ctype.from.LL = ifelse(ctype.from.LL == "GMP", "CMP", ctype.from.LL))
  # table(dat.meta.newannot.reannotate$ctype.from.LL)
# 
  # # dat.meta.reannotate <- subset(dat.meta.reannotate, select = -ctype.from.LL.secondbest)
  # dat.meta.newannot.reannotate$colcode <- sapply(dat.meta.newannot.reannotate$ctype.from.LL, function(x) colors.hash[[x]])
  # 
  # cell2ctypeLL <- hash::hash(dat.meta.newannot.reannotate$cell, dat.meta.newannot.reannotate$ctype.from.LL)
  
  dat.meta.reannotate <- dat.meta.final %>%
    rowwise() %>%
    mutate(ctype.from.LL = AssignHash(cell, closest.cells.hash, ctype), 
           ctype.from.LL = ifelse(ctype.from.LL == "CMP/GMP", "CMP", ctype.from.LL), 
           ctype.from.LL = ifelse(ctype.from.LL == "LT", "HSCs", ctype.from.LL), 
           ctype.from.LL = ifelse(ctype.from.LL == "ST", "HSCs", ctype.from.LL),
           colcode = AssignHash(ctype.from.LL, colors.hash, NA))
  table(dat.meta.reannotate$ctype.from.LL)
  
  # add to existing UMAP 
  m1 <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New + Old")) + 
    facet_wrap(~ctype.from.LL) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1)
  
  m2 <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New + Old")) + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2)
  
  m1.new <- ggplot(dat.meta.reannotate %>% filter(batch == "New"), aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New only")) + 
    facet_wrap(~ctype.from.LL) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.new)
  
  m2.new <- ggplot(dat.meta.reannotate %>% filter(batch == "New"), aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New only")) + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2.new)
  
  m1.old <- ggplot(dat.meta.reannotate %>% filter(batch == "Old"), aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "Old only")) + 
    facet_wrap(~ctype.from.LL) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.old)
  
  m2.old <- ggplot(dat.meta.reannotate %>% filter(batch == "Old"), aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "Old only")) + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2.old)
  
  # check DCs
  jctypes <- sort(unique(dat.meta.reannotate$ctype.from.LL))
  jbatches <- c("New", "Old")
  for (jctype in jctypes){
    print(jctype)
    for (jbatch in jbatches){
      print(jbatch)
      
      jsub <- dat.meta.reannotate %>% 
        mutate(highlight.ctype = ctype == jctype & batch == jbatch,
               highlight.ctype.LL = ctype.from.LL == jctype & batch == jbatch)
      m1.new2 <- ggplot(jsub %>%
                          arrange(highlight.ctype, highlight.ctype.LL), 
                        aes(x = umap1, y = umap2, color = highlight.ctype)) +
        geom_point() + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype, jbatch)) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      m1.new2.LL <- ggplot(jsub %>%
                             arrange(highlight.ctype.LL), 
                           aes(x = umap1, y = umap2, color = highlight.ctype.LL)) +
        geom_point() + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype, jbatch)) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      
      # ggplot(jsub %>% filter(highlight.ctype.LL), aes(x = umap1, y = umap2, color = ctype)) + 
      #   geom_point() + 
      #   facet_wrap(~ctype) + 
      #   theme_bw() + 
      #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      # 
      # ggplot(jsub %>% filter(highlight.ctype.LL), aes(x = umap1, y = umap2, color = ctype.from.LL.secondbest)) + 
      #   geom_point() + 
      #   facet_wrap(~ctype.from.LL.secondbest) + 
      #   theme_bw() + 
      #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      JFuncs::multiplot(m1.new2, m1.new2.LL, cols = 2)
      
    }
  }

  
   

  # Add colors --------------------------------------------------------------
  
  
  m1.final <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New + Old")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.final)
  
  m1.final.batch <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    facet_wrap(~batch) + 
    ggtitle(paste(jmark, "New + Old")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.final.batch)
  
  
  # Write outputs -----------------------------------------------------------
  
  fwrite(dat.meta.reannotate, file = outmeta, sep = "\t")
  dev.off()
}



# 
# # annotate some of the celltypes
# 
# dat.meta.to.change <- subset(dat.meta.reannotate, ctype %in% c("AllCells", "IL7RLinNeg", "LinNeg", "LSK", "HSPCs"))  %>%
#   dplyr::mutate(ctype = ctype.from.LL)
# dat.meta.keep <- subset(dat.meta.reannotate, !ctype %in% c("AllCells", "IL7RLinNeg", "LinNeg", "LSK", "HSPCs"))
# 
# dat.meta.changed <- rbind(dat.meta.to.change, dat.meta.keep)
# 
# m2 <- ggplot(dat.meta.changed, aes(x = umap1, y = umap2, color = ctype)) +
#   geom_point() + 
#   theme_bw() + 
#   facet_wrap(~ctype) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# JFuncs::multiplot(m1, m2, cols = 2)







