# Jake Yeung
# Date of Creation: 2022-05-15
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/35-UMAP_fix_colors.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks


# Load meta ----------------------------------------------------------------

# Old DC: #F0E442
# Old monocytes: #d8ce40

# New DC: #916559
# New monocytes: #F0E442

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt") 
  dat.meta <- fread(inf.meta) %>%
    rowwise() %>%
    mutate(colcodenew = ifelse(ctype.from.LL == "DCs", "#AF7C6E", colcode),
           colcodenew = ifelse(ctype.from.LL == "Monocytes", "#F0E422", colcodenew))
    # mutate(colcode = gsub("#d8ce40", "#41D676", colcode))
    # mutate(colcode = gsub("#d8ce40", "#82D641", colcode))
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcodenew))
dat.meta.colors.add <- data.frame(ctype.from.LL = c("GMP", "Tcells"), colcodenew = c("#d88543", "#04A804"))
dat.meta.colors <- rbind(dat.meta.colors, dat.meta.colors.add)

dat.meta.colors.unique <- dat.meta.colors[!duplicated(dat.meta.colors), ]

ctype2col <- hash::hash(dat.meta.colors.unique$ctype.from.LL, dat.meta.colors.unique$colcodenew)



# clst2color.old["GMP"] <- "#d88543"  # same as Granus
# clst2color.old["Tcells"] <- "#04A804"


# Show some ctypes --------------------------------------------------------


# ctypes.show <- c("Granulocytes", "DCs", "Monocytes", "pDCs")
ctypes.show <- list("Progs" = c("HSCs", "LT", "ST", "MPPs", "CMP", "MEP", "GMP"),
                    "Myeloids" = c("Granulocytes", "Monocytes", "DCs"), 
                    "Lymphoids" = c("Bcells", "NKs"), 
                    "Tcells" = c("Tcells"),
                    "pDCs" = c("pDCs"),
                    "Eryths" = c("MEP", "Eryths"))


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed"
dir.create(outdir)

for (jmark in jmarks){
  # ctypes.show <- sort(unique(subset(dat.meta.lst[[jmark]], batch == "New")$ctype.from.LL))
  # names(ctypes.show) <- ctypes.show
  
  dat.meta.sub <- dat.meta.lst[[jmark]] %>%
    rowwise() %>%
    mutate(ctype = ifelse(ctype == "GMP", "CMP", ctype))
  
  m <- ggplot(dat.meta.sub, aes(x = umap1, y = umap2, color = colcodenew)) +
    geom_point(alpha = 1) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcodenew,
                          guide = "legend") + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  outpdf <- file.path(outdir, paste0("ctypes_on_umap_groups_new_colors.", jmark, ".", Sys.Date(), ".pdf"))
  
  pdf(outpdf, useDingbats = FALSE)
  print(m)
  
  # show a Neutrophil just for example 
  (jcell <- subset(dat.meta.sub, ctype.from.LL == "Granulocytes")$cell[[1]])
  dat.meta.sub.example <- dat.meta.sub %>%
    rowwise() %>%
    mutate(colcodenew = ifelse(cell == jcell, colcodenew, "grey75"),
           is.cell = cell == jcell) %>%
    arrange(desc(colcodenew))
  
  m.ex <- ggplot(dat.meta.sub.example, aes(x = umap1, y = umap2, color = colcodenew, size = is.cell)) +
    geom_point(alpha = 1) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcodenew,
                          guide = "legend") + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.ex)
  
  print(m)
  for(jctype.name in names(ctypes.show)){
    
    jctype.vec <- ctypes.show[[jctype.name]]
    print(jctype.name)
    
    jtitle <- paste(jmark, jctype.vec)
    
    dat.meta.sub.label <- dat.meta.sub %>%
      rowwise() %>% 
      mutate(highlight.ctype = ctype %in% jctype.vec & batch == "New",
             col.highlight = ifelse(highlight.ctype, colcodenew, "gray75"), 
             col.byctype = ifelse(highlight.ctype, ctype2col[[ctype]], "gray75"),
             jalpha = ifelse(highlight.ctype, 1, 0.5)) %>%
      arrange(highlight.ctype)
    
    
    if (all(!dat.meta.sub.label$highlight.ctype)){
      print(paste("Nothing labeled for", jctype.vec, "skipping"))
      next
    }
    
    m2 <- ggplot(dat.meta.sub.label, aes(x = umap1, y = umap2, color = col.highlight, alpha = jalpha)) +
      geom_point() + 
      theme_bw() + 
      ggtitle(jctype.name) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcodenew,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m2)
    
    m2.batch <- ggplot(dat.meta.sub.label, aes(x = umap1, y = umap2, color = col.highlight, alpha = jalpha)) +
      geom_point() + 
      theme_bw() + 
      facet_wrap(~batch) + 
      ggtitle(jctype.name) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcodenew,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m2.batch)
    
    m3 <- ggplot(dat.meta.sub.label, aes(x = umap1, y = umap2, color = col.byctype, alpha = jalpha)) +
      geom_point() + 
      theme_bw() + 
      ggtitle(jctype.name) + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcodenew,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m3)
  }
  dev.off()
}


# Save outputs ------------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  fwrite(dat.meta.lst[[jmark]], file = file.path(outdir, paste0("umap_metadata_color_DC_monocyte_fixed.", jmark, ".", Sys.Date(), ".txt")), sep = "\t")
}

fwrite(dat.meta.colors.unique, file = file.path(outdir, paste0("dat_colors_DC_monocyte_fixed.", Sys.Date(), ".txt")), sep = "\t")

